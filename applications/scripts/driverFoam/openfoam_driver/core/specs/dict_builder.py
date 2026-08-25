#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     dict_builder
#
# Description
#     Constructs OpenFOAM dictionaries from programmatic specifications.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Solver-neutral primitives for synthesising OpenFOAM dictionary text.

This module owns the parts of dictionary synthesis that do not name any
particular solver's vocabulary:

- entry selection against a caller-supplied pool (`select_applicable_entries`),
- required-field checking (`check_required`),
- value resolution with `typical_value` fallback (`populate_values`),
- nested OpenFOAM block emission (`_set_nested`, `_serialize_block`) and
  value tokenisation (`_openfoam_value_token`),
- catalog-membership checks for caller-supplied override paths
  (`is_known_override_driver_path`).

The solver-specific builders that compose these primitives into concrete
`constant/` dictionaries — together with whatever scope-sentinel convention
their catalog uses — belong to the plugin that owns that vocabulary; for
cardiacFoam they live in `openfoam_driver.plugins.cardiacfoam.dict_builder`.
"""
from __future__ import annotations

from typing import Any

from openfoam_driver.dict_entries import DictEntry
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.core.specs.validation import (
    _entry_is_applicable,
    _predicate_matches,
    primary_phase,
    slot_key,
)


def select_applicable_entries(
    context: dict[str, Any],
    *,
    entries: list[DictEntry],
) -> list[DictEntry]:
    """Return only entries whose `applicable_when` predicate matches the
    context and whose `forbidden_when` predicate does NOT match. Entries
    with no `applicable_when` are always included. The `entries` pool is
    always supplied by the caller — this module knows no solver's catalog.
    Plugins may wrap this with their own default pool."""
    return [
        e for e in entries
        if _entry_is_applicable(e, context)
        and not any(
            _predicate_matches(context, key, expected)
            for key, expected in e.forbidden_when.items()
        )
    ]


def _is_required_in_context(
    entry: DictEntry,
    context: dict[str, Any],
) -> bool:
    """Decide whether an entry is required *for this context*.

    Three cases:
    - `required=True` AND `required_when` empty → always required.
    - `required=True` AND `required_when` non-empty → required ONLY when at
      least one `required_when` predicate matches. This reads the two fields
      together as the entry author's intent ("required, but only under
      these conditions"). Without this rule, an entry that is required only
      for one solver variant would fire missing-required errors on every
      other variant too.
    - `required=False` AND `required_when` non-empty → required only when a
      predicate matches (the validator's existing semantics).
    """
    if entry.required_when:
        return any(
            _predicate_matches(context, key, expected)
            for key, expected in entry.required_when.items()
        )
    return entry.required


def check_required(
    entries: list[DictEntry],
    populated: dict[str, str],
    *,
    context: dict[str, Any] | None = None,
) -> None:
    """Raise `ValueError` if any required entry in `entries` is missing from
    `populated`. Inapplicable entries are assumed already filtered out by
    `select_applicable_entries`; optional entries are silently ignored.

    `dynamic_path=True` entries are skipped — they describe template paths
    (e.g. ``domainCouplings.<name>.electroDomainCoupler``) rather than
    concrete required leaves, and this generic, plugin-agnostic function has
    no way to discover which concrete ``<name>`` instances a given run
    actually configures. Required-field enforcement for those concrete
    instances, if any, is a plugin concern: see e.g. the cardiacfoam
    plugin's ``_evaluate_dynamic_required_fields``
    (``plugins/cardiacfoam/validation.py``), which is not guaranteed to
    exist for every plugin's dynamic-path entries.

    The optional `context` enables `_is_required_in_context` to honour
    `required_when` predicates; when omitted, `required=True` is treated
    unconditionally for backward compat with simple callers.
    """
    missing: list[str] = []
    ctx = context if context is not None else {}
    for entry in entries:
        if entry.dynamic_path:
            continue
        if not _is_required_in_context(entry, ctx):
            continue
        key = slot_key(entry.driver_path)
        if key not in populated or populated[key] in (None, ""):
            missing.append(entry.driver_path)
    if missing:
        raise ValueError(
            "check_required: required entries have no value and "
            "no typical_value fallback was applicable:\n  - "
            + "\n  - ".join(missing)
        )


def populate_values(
    entries: list[DictEntry],
    context: dict[str, Any],
    *,
    typical_value_fallback: bool = True,
) -> dict[str, str]:
    """For each entry, resolve the final value to write into the dict.

    Precedence per entry:
      1. Explicit value already in `context` (from selectors or overrides).
      2. `entry.typical_value` if non-empty AND `typical_value_fallback`.
      3. Omit — caller's downstream `check_required` decides whether that's
         a problem for required entries.
    """
    import re
    populated: dict[str, str] = {}

    dynamic_entries = []
    for entry in entries:
        if getattr(entry, "dynamic_path", False):
            template = slot_key(entry.driver_path)
            # ANY <placeholder> is a wildcard, not just <name>/<electrode>.
            # Hardcoding those two silently dropped every override whose
            # template used a different placeholder -- including the ionic
            # constant overrides, whose <AC_name> segment never matched, so a
            # driver-written drug/channelopathy override emitted nothing at
            # all while validation reported success.
            pattern = _PLACEHOLDER_RE.sub(r"([^.]+)", re.escape(template))
            dynamic_entries.append((entry, template, re.compile(f"^{pattern}$")))

    active_instances: dict[str, set[tuple[str, ...]]] = {}
    for key, val in context.items():
        if val in (None, ""):
            continue
        for entry, template, regex in dynamic_entries:
            match = regex.match(key)
            if match:
                groups = match.groups()
                prefix = template.split(".<")[0]
                active_instances.setdefault(prefix, set()).add(groups)

    for entry in entries:
        if getattr(entry, "dynamic_path", False):
            template = slot_key(entry.driver_path)
            prefix = template.split(".<")[0]
            if prefix in active_instances:
                for groups in active_instances[prefix]:
                    # Substitute captured groups positionally, so a template
                    # with any number of placeholders reconstructs correctly.
                    concrete_key = template
                    for captured in groups:
                        concrete_key = _PLACEHOLDER_RE.sub(captured, concrete_key, count=1)

                    if concrete_key in context and context[concrete_key] not in (None, ""):
                        populated[concrete_key] = str(context[concrete_key])
                    elif typical_value_fallback and entry.typical_value:
                        populated[concrete_key] = entry.typical_value
            continue

        key = slot_key(entry.driver_path)
        if key in context and context[key] not in (None, ""):
            populated[key] = str(context[key])
            continue
        if typical_value_fallback and entry.typical_value:
            conflict = False
            for mx_path in getattr(entry, "mutually_exclusive_with", ()):
                mx_key = slot_key(mx_path)
                if mx_key in context and context[mx_key] not in (None, ""):
                    conflict = True
                    break
            if not conflict:
                populated[key] = entry.typical_value
            continue
    return populated


def _populated_to_run(
    populated: dict[str, str],
    entries: list[DictEntry],
) -> RunDocument:
    """Distribute populated values into a Run document keyed by each
    entry's primary phase. Selector keys (which may not correspond to any
    entry, but always do here for the dict-builder entry pool) are placed
    in the physics slice as a sensible default."""
    config: dict[str, dict[str, str]] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    placed: set[str] = set()
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key not in populated:
            continue
        ph = primary_phase(entry) or "physics"
        config[ph][key] = populated[key]
        placed.add(key)
    # Any populated keys without a matching entry land in physics. This
    # only triggers for selector keys that don't correspond to DictEntry —
    # uncommon, but safe.
    for key, val in populated.items():
        if key not in placed:
            config["physics"][key] = val
    return RunDocument(id="dict_builder", name="dict_builder",
                       status="draft", config=config)


import re as _re

_PLACEHOLDER_RE = _re.compile(r"<[A-Za-z_][A-Za-z0-9_]*>")


def is_known_override_driver_path(
    key: str,
    *,
    driver_context=None,
) -> bool:
    """True if `key` matches a real dict-entry driver_path in the active
    plugin's catalog.

    Matching is prefix-agnostic (via `slot_key`) and honours `dynamic_path`
    templates (e.g. ``domainCouplings.<name>.electroDomainCoupler``) by
    treating any ``<placeholder>`` segment as a wildcard. This is a pure
    catalog-membership check — it does not consider whether the entry is
    *applicable* in a given selector context (that's `select_applicable_entries`'s
    job, run later inside the plugin's builder). Used by callers that
    need to reject a caller-supplied override path outright before it is ever
    passed to `build_and_launch` (e.g. `sweep_routing.route_case_values`),
    rather than silently accepting an override that has no matching entry
    anywhere and therefore no effect.
    """
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    normalized = slot_key(key)
    for entry in driver_context.capabilities.dictionaries.catalog().entries:
        entry_key = slot_key(entry.driver_path)
        if getattr(entry, "dynamic_path", False):
            # re.escape leaves `<`, `>`, letters and `_` untouched, so
            # escaping first and substituting placeholders after is safe.
            pattern = _PLACEHOLDER_RE.sub(r"[^.]+", _re.escape(entry_key))
            if _re.fullmatch(pattern, normalized):
                return True
        elif entry_key == normalized:
            return True
    return False


def _set_nested(node: dict, path: list[str], value: Any) -> None:
    """Insert `value` at `path` inside the nested dict `node`, creating
    intermediate sub-dicts as needed. A leaf already present is
    overwritten — the populated dict has unique slot_keys so this is safe."""
    cursor = node
    for segment in path[:-1]:
        cursor = cursor.setdefault(segment, {})
        # If a prior leaf collided with a sub-block name, replace the leaf
        # with a sub-block — should not happen with the current catalog but
        # is defensive.
        if not isinstance(cursor, dict):
            raise ValueError(
                f"Path collision in serialiser at segment {segment!r}: a leaf "
                f"value exists where a sub-block is needed."
            )
    cursor[path[-1]] = value


def _openfoam_value_token(value: str) -> str:
    """Return `value` in a form OpenFOAM's tokenizer accepts as a dict value.

    OpenFOAM lexes a bare token starting with a digit as a number. A word like
    ``3D`` therefore reads as label ``3`` followed by junk, and the run dies
    with "expected word, found label 3" -- which is exactly what made
    ``$ELECTRO_MODEL_COEFFS.dimension`` unusable through the driver. Tutorials
    write ``dimension "3D";``.

    Quote ONLY that case. Anything that parses as a number, is already quoted,
    or is a compound token (vector, list, dimension set) must pass through
    untouched -- quoting those would break dictionaries that work today.
    """
    if not isinstance(value, str) or not value:
        return value
    token = value.strip()
    if not token or not token[0].isdigit():
        return value
    if any(ch.isspace() for ch in token) or token[0] in "([{\"":
        return value
    try:
        float(token)
    except ValueError:
        return f'"{token}"'
    return value


def _serialize_block(tree: dict, indent: int) -> str:
    """Emit nested OpenFOAM block syntax. Leaves are `key value;`,
    sub-blocks are `key\\n{\\n  ...\\n}` recursively."""
    lines: list[str] = []
    pad = " " * indent
    for key, value in tree.items():
        if isinstance(value, dict):
            lines.append(f"{pad}{key}")
            lines.append(f"{pad}{{")
            lines.append(_serialize_block(value, indent + 4))
            lines.append(f"{pad}}}")
        else:
            lines.append(f"{pad}{key} {_openfoam_value_token(value)};")
    return "\n".join(lines)
