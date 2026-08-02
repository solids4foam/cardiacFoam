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

"""Scratch-construct OpenFOAM dictionary files from agent intent.

`build_electro_properties` synthesises a complete `constant/electroProperties`
text from selectors + overrides. The pipeline reuses the existing dict-entry
catalog (`dict_entries.py`), the structured-constraint validator
(`validation.py`), and the path conventions encoded in `slot_key`.

The builder composes existing primitives. Every output passes through `validate_run`
before being returned; an agent that gets a string back is guaranteed it
is validator-clean.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from openfoam_driver.dict_entries import (
    DictEntry,
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import (
    _entry_is_applicable,
    _predicate_matches,
    primary_phase,
    slot_key,
    validate_run,
)


def _all_electro_entries() -> list[DictEntry]:
    ordered_keys = [
        "top_level",
        "common_model_coeffs",
        "monodomain",
        "bidomain",
        "bath_potential_domain",
        "eikonal_diffusion",
        "ionic_heterogeneity",
        "ionic_constant_overrides",
        "batched_integrator",
        "active_tension",
        "ode_solver_passthrough",
        "single_cell_stimulus",
        "conduction_system",
        "domain_couplings",
        "ecg"
    ]
    out: list[DictEntry] = []
    for k in ordered_keys:
        if k in get_electro_property_entry_groups():
            out.extend(get_electro_property_entry_groups()[k])
    for k, group in get_electro_property_entry_groups().items():
        if k not in ordered_keys:
            out.extend(group)
    return out


# Block-presence virtual keys auto-inferred from override prefixes.
# When an agent declares any override under `bathPotentialDomain.*`,
# `ecgDomains.*`, or `conductionNetworkDomains.*`, the dict_builder sets
# the corresponding `$..._present`/`$..._configured` virtual key on the
# resolved context. Entries gated by `applicable_when={"$..._...": True}`
# then become visible and their typical_value fallbacks fire. Without an
# override under one of these prefixes the block stays off — plain
# bidomain does not get bath leaves leaking in.
_VIRTUAL_PRESENCE_TRIGGERS: tuple[tuple[str, str], ...] = (
    ("bathPotentialDomain.", "$bathPotentialDomain_configured"),
    ("ecgDomains.", "$ecgDomains_present"),
    ("conductionNetworkDomains.", "$conductionNetworkDomains_present"),
)


def _infer_virtual_presence(ctx: dict[str, Any]) -> None:
    """Set block-presence virtual keys in-place on `ctx` whenever any
    real slot_key starts with one of the documented block prefixes.
    Idempotent and safe to call after the override merge."""
    for prefix, virtual_key in _VIRTUAL_PRESENCE_TRIGGERS:
        if virtual_key in ctx:
            continue
        for existing_key in ctx:
            if existing_key.startswith(prefix):
                ctx[virtual_key] = True
                break

    # Support OR logic for ionicHeterogeneity applicability
    from openfoam_driver.dict_entries import get_heterogeneity_models
    if ctx.get("myocardiumSolver") == "eikonalSolver" or ctx.get("ionicModel") in get_heterogeneity_models():
        ctx["$ionicHeterogeneity_supported"] = True


def resolve_context(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Collapse selectors + overrides into a single `{slot_key: value}` dict.

    Selectors enter at their raw key (`myocardiumSolver`, `ionicModel`, ...);
    overrides go through `slot_key` so the `$ELECTRO_MODEL_COEFFS.` prefix
    is stripped. This is the same context shape `validation._flatten_context`
    produces from a Run document — so the validator can be reused unchanged.

    Also infers virtual presence keys (`$bathPotentialDomain_configured`,
    `$ecgDomains_present`, `$conductionNetworkDomains_present`) so the
    matching `applicable_when` predicates fire only when the agent has
    actually declared overrides under the corresponding block.
    """
    ctx: dict[str, Any] = dict(selectors)
    if overrides:
        for driver_path, value in overrides.items():
            ctx[slot_key(driver_path)] = value
    if ctx.get("myocardiumSolver") in {
        "monodomainSolver", "eikonalSolver", "bidomainSolver"
    }:
        ctx.setdefault("conductivitySource", "uniform")
    _infer_virtual_presence(ctx)
    return ctx


def select_applicable_entries(
    context: dict[str, Any],
    *,
    entries: list[DictEntry] | None = None,
) -> list[DictEntry]:
    """Return only entries whose `applicable_when` predicate matches the
    context and whose `forbidden_when` predicate does NOT match. Entries
    with no `applicable_when` are always included. The `entries` kwarg lets
    callers (and tests) curate the input pool; default is every electro entry."""
    pool = entries if entries is not None else _all_electro_entries()
    return [
        e for e in pool
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
      these conditions"). Without this rule, entries like `phiERefPoint`
      (required=True, required_when={"myocardiumSolver":"bidomainSolver"})
      would fire missing-required errors on singleCell cases.
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
    concrete required leaves. The user's overrides supply concrete paths
    when those blocks are actually configured; the cross-block validators catch
    dangling references at run-construction time.

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
            "build_electro_properties: required entries have no value and "
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
            pattern = re.escape(template).replace("<name>", r"([^\.]+)")
            pattern = pattern.replace("<electrode>", r"([^\.]+)")
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
                    concrete_key = template
                    if "<name>" in concrete_key and len(groups) > 0:
                        concrete_key = concrete_key.replace("<name>", groups[0])
                    if "<electrode>" in concrete_key and len(groups) > 1:
                        concrete_key = concrete_key.replace("<electrode>", groups[-1])
                    elif "<electrode>" in concrete_key and len(groups) == 1:
                        concrete_key = concrete_key.replace("<electrode>", groups[0])

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


def _foamfile_preamble(object_name: str) -> str:
    """Standard FoamFile preamble for a given `object` value.

    The object name distinguishes the two dict files this module writes —
    `electroProperties` vs `physicsProperties`. Future dict targets follow
    the same pattern.
    """
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | cardiacFoam dict_builder synthesis              |\n"
        "|   \\\\  /    A nd           |                                                 |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        "    location    \"constant\";\n"
        f"    object      {object_name};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


# Backwards-compat alias — some callers and tests reference this name.
_FOAMFILE_PREAMBLE = _foamfile_preamble("electroProperties")


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


def _check_no_forbidden_selectors(context: dict[str, Any]) -> None:
    """Raise ValueError if the caller explicitly set a key that is forbidden
    in the current context.  'Explicitly set' means the slot key appears in
    *context* — i.e., the caller passed it as a selector or override."""
    for entry in _all_electro_entries():
        if not entry.forbidden_when:
            continue
        key = slot_key(entry.driver_path)
        if key not in context:
            continue
        for pred_key, expected in entry.forbidden_when.items():
            if _predicate_matches(context, pred_key, expected):
                raise ValueError(
                    f"build_electro_properties: '{key}' is forbidden when "
                    f"{pred_key}={context.get(pred_key)!r}."
                )


def build_electro_properties(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
    typical_value_fallback: bool = True,
) -> str:
    """Synthesise a complete `electroProperties` dict from intent.

    Args:
        selectors: top-level discriminators (myocardiumSolver, ionicModel,
            tissue, ...). Required keys depend on the chosen solver.
        overrides: full driver_path → value mappings for entries whose
            `typical_value` is not appropriate.
        typical_value_fallback: when True (default), applicable entries
            with no override fall back to `DictEntry.typical_value` if
            declared. When False, only explicit overrides count.

    Returns:
        OpenFOAM-format text including the standard `FoamFile` preamble.

    Raises:
        ValueError: required+applicable entry has no value, mutex violation,
            or any structured-constraint violation from `validate_run`.
    """
    context = resolve_context(selectors, overrides=overrides)

    # Pre-check: if the caller explicitly provides a key that is forbidden
    # in this context, reject immediately rather than silently ignoring it.
    _check_no_forbidden_selectors(context)

    entries = select_applicable_entries(context)
    populated = populate_values(
        entries, context, typical_value_fallback=typical_value_fallback,
    )

    # Safety net: run the full validator scoped to electro entries only.
    # The validator now subsumes required-field checks (its section 1
    # honours `required_when` + `dynamic_path` the same way `check_required`
    # does), so we don't pre-call `check_required` from the public builder
    # entry-point. `check_required` stays exported for callers that want
    # just the required-field subset.
    run = _populated_to_run(populated, entries)
    errors = [e for e in validate_run(run, entries=entries) if e.level == "error"]
    if errors:
        raise ValueError(
            "build_electro_properties: validator rejected synthesised dict:\n  - "
            + "\n  - ".join(e.message for e in errors)
        )

    body = _serialize(populated, entries, selectors["myocardiumSolver"])
    return _FOAMFILE_PREAMBLE + "\n" + body


# Slot keys that map to selectors (top-level discriminators), not overrides.
_SELECTOR_KEYS: frozenset[str] = frozenset(
    {"myocardiumSolver", "ionicModel", "tissue", "conductivitySource"}
)
SELECTOR_KEYS: frozenset[str] = _SELECTOR_KEYS  # public alias for external consumers (e.g. sweep_routing.py)

import re as _re

_PLACEHOLDER_RE = _re.compile(r"<[a-zA-Z_]+>")


def is_known_override_driver_path(key: str) -> bool:
    """True if `key` matches a real electro or physics dict-entry driver_path.

    Matching is prefix-agnostic (via `slot_key`) and honours `dynamic_path`
    templates (e.g. ``domainCouplings.<name>.electroDomainCoupler``) by
    treating any ``<placeholder>`` segment as a wildcard. This is a pure
    catalog-membership check — it does not consider whether the entry is
    *applicable* in a given selector context (that's `select_applicable_entries`'s
    job, run later inside `build_electro_properties`). Used by callers that
    need to reject a caller-supplied override path outright before it is ever
    passed to `build_and_launch` (e.g. `sweep_routing.route_case_values`),
    rather than silently accepting an override that has no matching entry
    anywhere and therefore no effect.
    """
    normalized = slot_key(key)
    for entry in list(_all_electro_entries()) + list(PHYSICS_PROPERTY_ENTRIES):
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

_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."


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
            lines.append(f"{pad}{key} {value};")
    return "\n".join(lines)


def _serialize(
    populated: dict[str, str],
    entries: list[DictEntry],
    myocardium_solver: str,
) -> str:
    """Group populated values by scope and emit the OpenFOAM dict body.

    Top-level keys (entries whose `driver_path` does not start with the
    `$ELECTRO_MODEL_COEFFS.` prefix) are emitted at the root. Everything
    else nests under the resolved `<solver>Coeffs` block.
    """
    import re
    top_level: dict[str, str] = {}
    coeffs: dict = {}

    dynamic_patterns = []
    for entry in entries:
        if getattr(entry, "dynamic_path", False):
            template = slot_key(entry.driver_path)
            pattern = re.escape(template).replace("<name>", r"([^\.]+)")
            pattern = pattern.replace("<electrode>", r"([^\.]+)")
            dynamic_patterns.append((entry, re.compile(f"^{pattern}$")))

    for concrete_key, value in populated.items():
        comment = ""
        matched_entry = None
        for entry in entries:
            if not getattr(entry, "dynamic_path", False):
                if slot_key(entry.driver_path) == concrete_key:
                    matched_entry = entry
                    break

        if not matched_entry:
            for entry, regex in dynamic_patterns:
                if regex.match(concrete_key):
                    matched_entry = entry
                    break

        if matched_entry and matched_entry.driver_path.startswith(_COEFFS_PREFIX):
            segments = concrete_key.split(".")
            _set_nested(coeffs, segments, value)
        else:
            top_level[concrete_key] = value

    parts: list[str] = []
    for key, value in top_level.items():
        parts.append(f"{key} {value};")
    if coeffs:
        coeffs_scope = f"{myocardium_solver}Coeffs"
        parts.append("")
        parts.append(coeffs_scope)
        parts.append("{")
        parts.append(_serialize_block(coeffs, indent=4))
        parts.append("}")
    parts.append("")
    return "\n".join(parts)


def _entry_scope_and_key(
    driver_path: str,
    coeffs_scope: str,
) -> tuple[list[str] | None, str]:
    """Resolve a driver_path to (scope_path, key) suitable for read_foam_entry.

    Examples:
        "myocardiumSolver"                          → (None, "myocardiumSolver")
        "$ELECTRO_MODEL_COEFFS.solutionAlgorithm"   → (["monodomainSolverCoeffs"], "solutionAlgorithm")
        "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude"
            → (["singleCellSolverCoeffs", "singleCellStimulus"], "stim_amplitude")
    """
    if driver_path.startswith(_COEFFS_PREFIX):
        segments = driver_path[len(_COEFFS_PREFIX):].split(".")
        if len(segments) == 1:
            return [coeffs_scope], segments[0]
        return [coeffs_scope] + segments[:-1], segments[-1]
    return None, driver_path


def parse_electro_properties(
    electro_properties_path: "Any",
) -> dict[str, dict[str, str]]:
    """Parse an existing ``electroProperties`` file into selectors + overrides.

    Reads the file using the same ``get_electro_property_entry_groups()`` catalog
    that :func:`build_electro_properties` writes from. Returns a dict that
    round-trips through :func:`build_electro_properties`.

    - ``selectors``: ``myocardiumSolver``, ``ionicModel``, ``tissue`` (when
      present and applicable).
    - ``overrides``: full ``driver_path → value`` for every non-default entry
      found in the file. Values equal to ``entry.typical_value`` are omitted
      (the builder fills them automatically). Entries with ``dynamic_path=True``
      are skipped entirely.
    - ``ignored_keys``: the ``driver_path`` of every catalog entry this parser
      structurally does not round-trip (the ``dynamic_path=True`` families). If
      the source dict sets any of these, they will NOT reappear on a rebuild —
      inspect the source dict manually. Surfacing them here replaces the old
      silent drop. (Keys entirely outside the catalog remain un-enumerated —
      the parser only reads catalogued paths.)

    Returns:
        ``{"selectors": {...}, "overrides": {...}, "ignored_keys": [...]}``
    """
    from pathlib import Path as _Path
    from openfoam_driver.specs.common import detect_myocardium_solver_name
    from openfoam_driver.core.runtime.mutators import read_foam_entry

    electro_properties_path = _Path(electro_properties_path)
    solver = detect_myocardium_solver_name(electro_properties_path)
    coeffs_scope = f"{solver}Coeffs"

    selectors: dict[str, str] = {"myocardiumSolver": solver}
    overrides: dict[str, str] = {}
    ignored_keys: list[str] = []

    for entry in _all_electro_entries():
        if entry.dynamic_path:
            ignored_keys.append(entry.driver_path)
            continue
        scope_path, key = _entry_scope_and_key(entry.driver_path, coeffs_scope)
        value = read_foam_entry(electro_properties_path, key, scope=scope_path)
        if value is None:
            continue
        sk = slot_key(entry.driver_path)
        if sk in _SELECTOR_KEYS:
            selectors[sk] = value
        elif value != entry.typical_value:
            overrides[entry.driver_path] = value

    return {
        "selectors": selectors,
        "overrides": overrides,
        "ignored_keys": sorted(set(ignored_keys)),
    }


def build_physics_properties(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
    typical_value_fallback: bool = True,
) -> str:
    """Synthesise a complete `physicsProperties` dict from intent.

    Mirrors :func:`build_electro_properties` but against
    :data:`PHYSICS_PROPERTY_ENTRIES`. There is no ``<solver>Coeffs``
    wrapper — every physics key lives at the dict root. Today the only
    entry is ``type``; future physics-level selectors slot in unchanged.

    Args:
        selectors: top-level physics keys (e.g. ``{"type": "electroModel"}``).
        overrides: full driver_path → value mappings for future expansion.
        typical_value_fallback: when True, applicable entries with no
            override fall back to ``DictEntry.typical_value`` if declared.

    Returns:
        OpenFOAM-format text with a ``physicsProperties``-typed FoamFile
        preamble.

    Raises:
        ValueError: required entry has no value, or any structured
            constraint violation from `validate_run`.
    """
    context = resolve_context(selectors, overrides=overrides)
    # Scope to physics entries — electro entries don't belong here.
    entries = select_applicable_entries(context, entries=list(PHYSICS_PROPERTY_ENTRIES))
    populated = populate_values(
        entries, context, typical_value_fallback=typical_value_fallback,
    )

    run = _populated_to_run(populated, entries)
    errors = [e for e in validate_run(run, entries=entries) if e.level == "error"]
    if errors:
        raise ValueError(
            "build_physics_properties: validator rejected synthesised dict:\n  - "
            + "\n  - ".join(e.message for e in errors)
        )

    # Physics keys are root-level — no <solver>Coeffs wrapper.
    body_lines: list[str] = []
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key in populated:
            body_lines.append(f"{key} {populated[key]};")
    body = "\n".join(body_lines) + "\n"
    return _foamfile_preamble("physicsProperties") + "\n" + body


def build_and_launch(
    electro_selectors: dict[str, str],
    *,
    physics_selectors: dict[str, str],
    case_dir: "Path",
    electro_overrides: "dict[str, str] | None" = None,
    physics_overrides: "dict[str, str] | None" = None,
    overwrite: bool = False,
    dry_run: bool = False,
    pre_solve_commands: "Sequence[str | Sequence[str]] | None" = None,
    openfoam_bashrc: "str | Path | None" = None,
    delta_t: "float | str | None" = None,
    end_time: "float | str | None" = None,
    dx: "float | None" = None,
) -> dict:
    """Build both dicts, write them to ``case_dir/constant/``, and (if
    not dry_run) launch the engine on the resulting case.

    Args:
        electro_selectors: selectors for build_electro_properties.
        physics_selectors: selectors for build_physics_properties.
        case_dir: target case directory. Will be created if absent.
        electro_overrides / physics_overrides: optional overrides per
            builder semantics.
        overwrite: when False (default), an existing
            ``case_dir/constant/electroProperties`` raises FileExistsError.
        dry_run: when True, writes the dicts and returns without
            running the engine.
        pre_solve_commands: optional commands to run in ``case_dir`` before
            ``cardiacFoam``. Each item is a string (shell-split) or a list of
            strings. Example: ``["vtkUnstructuredToFoam",
            "setTorsoOrganConductivityField"]``.
        openfoam_bashrc: when set, each command is run after sourcing this
            OpenFOAM bashrc.
        dx: mesh resolution (metres, isotropic cell size) for the generic
            default ``blockMeshDict`` provisioned for spatial solvers with no
            author-supplied mesh. Only meaningful for
            ``monodomainSolver``/``bidomainSolver``/``eikonalSolver``; raises
            ``ValueError`` for ``singleCellSolver`` (no spatial geometry)
            rather than silently having no effect, and if it does not evenly
            divide the default slab's fixed size (no silent rounding).
            Meaningless for real anatomical meshes imported via
            ``vtkUnstructuredToFoam`` -- this only controls the generic
            default slab.

    Returns:
        A dict carrying ``case_dir`` (str) and either
        ``status="dry_run_complete"`` or the engine result list under
        ``results``.

    Raises:
        FileExistsError: case_dir/constant/electroProperties exists and
            overwrite=False.
        ValueError: validator rejected one of the synthesized dicts.
    """
    from pathlib import Path as _Path

    case_dir = _Path(case_dir)
    constant_dir = case_dir / "constant"
    electro_path = constant_dir / "electroProperties"
    physics_path = constant_dir / "physicsProperties"

    if electro_path.exists() and not overwrite:
        raise FileExistsError(
            f"{electro_path} already exists; pass overwrite=True to replace."
        )

    electro_text = build_electro_properties(
        electro_selectors, overrides=electro_overrides,
    )
    physics_text = build_physics_properties(
        physics_selectors, overrides=physics_overrides,
    )

    constant_dir.mkdir(parents=True, exist_ok=True)
    electro_path.write_text(electro_text)
    physics_path.write_text(physics_text)

    system_dir = case_dir / "system"
    system_dir.mkdir(parents=True, exist_ok=True)

    from openfoam_driver.specs.system_templates import get_fv_schemes, get_fv_solution, build_control_dict
    myocardium_solver = electro_selectors.get("myocardiumSolver", "monodomainSolver")

    fv_schemes_path = system_dir / "fvSchemes"
    if not fv_schemes_path.exists() or overwrite:
        fv_schemes_path.write_text(get_fv_schemes(myocardium_solver))

    fv_solution_path = system_dir / "fvSolution"
    if not fv_solution_path.exists() or overwrite:
        fv_solution_path.write_text(get_fv_solution(myocardium_solver))

    control_dict_path = system_dir / "controlDict"
    if not control_dict_path.exists() or overwrite:
        dt = delta_t if delta_t is not None else 1e-4
        et = end_time if end_time is not None else 1.0
        control_dict_path.write_text(build_control_dict(delta_t=dt, end_time=et))

    if delta_t is not None or end_time is not None:
        from openfoam_driver.core.runtime.mutators import update_control_dict
        update_control_dict(
            case_dir / "system" / "controlDict",
            delta_t=delta_t,
            end_time=end_time,
        )

    from openfoam_driver.specs.mesh_provisioning import provision_mesh
    needs_block_mesh = provision_mesh(
        case_dir=case_dir, myocardium_solver=myocardium_solver, dx_m=dx,
    )

    if dry_run:
        return {
            "case_dir": str(case_dir),
            "status": "dry_run_complete",
            "needs_block_mesh": needs_block_mesh,
        }

    from openfoam_driver.specs.tutorials.generic_case import make_spec
    from openfoam_driver.core.runtime.execution_context import resolve_execution_context
    from openfoam_driver.core.runtime.openfoam_environment import load_openfoam_environment
    from openfoam_driver.core.runtime.workflow import normalize_workflow_dag, validate_workflow_commands
    from openfoam_driver.core.runtime.workflow_orchestrator import run_workflow
    from openfoam_driver.core.runtime.workflow_state import initial_workflow_state

    spec = make_spec(
        tutorials_root=case_dir.parent,
        case_dir_name=case_dir.name,
        solver_command="cardiacFoam",
        pre_solve_commands=list(pre_solve_commands or ()),
        openfoam_bashrc=openfoam_bashrc,
    )
    execution_context = resolve_execution_context(spec)
    workflow_dag, _dag_diagnostics = normalize_workflow_dag(spec.metadata.get("workflow_dag"))
    command_diagnostics = validate_workflow_commands(workflow_dag)
    if command_diagnostics:
        raise ValueError(
            "build_and_launch's workflow_dag failed command validation: "
            + "; ".join(d.message for d in command_diagnostics)
        )

    env = load_openfoam_environment(explicit_bashrc=openfoam_bashrc).env
    outcome = run_workflow(
        workflow_dag,
        initial_workflow_state(workflow_dag),
        case_root=execution_context.case_root,
        output_dir=execution_context.output_dir,
        env=env,
    )
    return {
        "case_dir": str(case_dir),
        "status": "complete" if outcome.state.status == "completed" else "failed",
        "workflow_state": outcome.state.to_json(),
    }
