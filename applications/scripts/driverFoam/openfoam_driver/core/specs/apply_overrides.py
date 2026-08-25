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
#     apply_overrides
#
# Description
#     Mechanically apply an agent-chosen override set to a case's dicts for
#     `step --strict --apply`. Validates each override for *applyability*
#     (catalog-addressable AND writable by the router) before any write, then
#     routes controlDict leaves to update_control_dict and "$TOKEN." leaves
#     through whichever OverrideScope the active plugin declares for that
#     token (see :class:`OverrideScope`). The driver never decides *what* to
#     change — the agent authors the override set; this only applies a
#     validated one.
#
#     Core knows no scope tokens itself: a plugin declares each one via
#     ``PluginCapabilities.override_scopes`` (the built-in cardiac plugin
#     declares exactly one, $ELECTRO_MODEL_COEFFS -> constant/electroProperties
#     -> plugins/cardiacfoam/overrides.py::electro_model_coeffs_scope). An
#     override whose token matches no declared scope is rejected with an
#     explicit "unknown scope token" error rather than the generic
#     "not catalog-addressable" this module used to raise for every $-prefixed
#     miss. specs/validation.py and scripts/_dict_keys_scanner.py separately
#     strip a leading "$TOKEN." for phase-slice/drift-scan normalization —
#     that's a syntactic transform needing no plugin lookup, so it stayed
#     generalized independently rather than routed through this capability.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path, PurePath
from typing import Any, Callable, Iterable

from ..runtime.mutators import update_foam_entry


@dataclass(frozen=True)
class OverrideScope:
    """One ``$TOKEN.`` override scope a plugin declares for `step --strict
    --apply`.

    token: the bare scope name after ``$`` and before the first ``.`` (e.g.
        ``"ELECTRO_MODEL_COEFFS"`` for ``"$ELECTRO_MODEL_COEFFS.myocardiumSolver"``).
    file_relpath: the case-relative dict file this scope's overrides write
        into (e.g. ``"constant/electroProperties"``).
    catalog_group: the dictionary-catalog group name this scope's overrides
        are validated against (``DictionaryCatalogCapability.catalog()
        .entries_for(catalog_group)``).
    resolve_entry: given the full ``driver_path`` and the case root, return
        ``(scope_path, key)`` ready for
        :func:`openfoam_driver.core.runtime.mutators.update_foam_entry`.
        Plugin-owned: how a token's dotted suffix maps onto nested OpenFOAM
        scope segments is catalog-specific (e.g. which ``<solver>Coeffs``
        block is active for this case), not something core can infer.
    """

    token: str
    file_relpath: str
    catalog_group: str
    resolve_entry: Callable[[str, Path], tuple[list[str] | None, str]]


@dataclass(frozen=True)
class RegenerationScope:
    """One bare (non-``$``-prefixed) "selector" override a plugin declares
    for `step --strict --apply` that must REGENERATE a dict file rather
    than key-patch it, because changing the value restructures the file --
    renames a sub-block, changes which sibling keys are legal -- instead of
    changing one leaf in place. The motivating case is cardiacFoam's
    ``myocardiumSolver``: switching it renames ``<oldSolver>Coeffs`` to
    ``<newSolver>Coeffs`` and flips which keys the catalog's
    ``applicable_when``/``required_when``/``forbidden_when`` predicates
    allow, none of which ``update_foam_entry`` (a single key/value/scope
    patch) can express.

    selector_keys: the bare ``driver_path`` names this scope owns (e.g.
        ``{"myocardiumSolver"}`` for the cardiac plugin). Deliberately a
        small, explicit set: a selector only belongs here if changing it
        actually restructures the file. Plain leaf selectors that only
        change a value in place (e.g. cardiacFoam's ``ionicModel``,
        ``tissue``) stay on the ordinary ``$TOKEN.`` :class:`OverrideScope`
        path instead -- routing them through regeneration too would be a
        capability the catalog does not need yet.
    file_relpath: the case-relative dict file this scope regenerates (e.g.
        ``"constant/electroProperties"``).
    catalog_group: the dictionary-catalog group name used to validate enum
        values for these selector keys (mirrors
        :attr:`OverrideScope.catalog_group`).
    regenerate: given the on-disk file path, the full ``driver_path``, the
        new value, and any OTHER ``$TOKEN.``-scoped overrides from the same
        `step --strict --apply` call that target this same
        ``file_relpath`` (``{driver_path: value}``, empty if none),
        rewrite the file in place from its current content with that one
        selector changed. The extra-overrides map exists because
        ``update_foam_entry`` can only patch a key that already exists --
        a solver switch can make a *new* key required with no catalog
        default (e.g. eikonalSolver's ``stimulusLocationMin``, absent from
        a monodomain source and un-defaultable, case-specific geometry),
        and there would otherwise be no way for such a value to reach the
        file: too late to key-patch it in afterward (nothing to patch),
        and the rebuild has no default to fall back on. Plugin-owned: how
        to decompose the existing file into selectors/overrides, rebuild
        it, and preserve whatever the rebuild pipeline cannot itself
        round-trip is catalog-specific, not something core can infer.
    """

    selector_keys: frozenset[str]
    file_relpath: str
    catalog_group: str
    regenerate: Callable[[Path, str, str, dict[str, str]], None]


def _is_safe_system_path(path_str: str) -> bool:
    """Validate that the path is strictly inside system/ and has no traversal segments."""
    if not path_str.startswith("system/"):
        return False
    path = PurePath(path_str)
    return not path.is_absolute() and ".." not in path.parts


class OverrideError(ValueError):
    """An override is malformed, non-applyable, out-of-enum, or failed to apply."""


def _catalog_entries(
    driver_context=None,
) -> tuple[set[str], dict[str, Any], tuple[OverrideScope, ...], tuple[RegenerationScope, ...]]:
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    catalog = driver_context.capabilities.dictionaries.catalog()
    scopes = driver_context.capabilities.override_scopes.scopes()
    regeneration_scopes = driver_context.capabilities.dict_regeneration.scopes()
    scoped_entries: dict[str, Any] = {}
    for scope in scopes:
        for entry in catalog.entries_for(scope.catalog_group):
            scoped_entries[entry.driver_path] = entry
    for regen_scope in regeneration_scopes:
        for entry in catalog.entries_for(regen_scope.catalog_group):
            scoped_entries.setdefault(entry.driver_path, entry)
    return (
        {entry.driver_path for entry in catalog.entries_for("controlDict")},
        scoped_entries,
        scopes,
        regeneration_scopes,
    )


def _match_dynamic_entry(dp: str, all_entries: Iterable[Any]) -> Any | None:
    """Return the dynamic catalog entry whose template matches concrete *dp*."""
    for entry in all_entries:
        if not getattr(entry, "dynamic_path", False):
            continue

        template = entry.driver_path
        pattern_parts: list[str] = []
        previous_end = 0
        for placeholder in re.finditer(r"<[^.<>]+>", template):
            pattern_parts.append(re.escape(template[previous_end:placeholder.start()]))
            pattern_parts.append(r"[^.]+")
            previous_end = placeholder.end()
        pattern_parts.append(re.escape(template[previous_end:]))

        if re.fullmatch("".join(pattern_parts), dp):
            return entry
    return None


def _scope_token(dp: str) -> str:
    """Return the bare token between "$" and the first "." (or the whole
    remainder if there is no "."). ``dp`` must already be known to start
    with "$"."""
    return dp[1:].split(".", 1)[0]


def validate_overrides(overrides: Any, *, driver_context=None) -> None:
    """Reject anything not safely applyable, *before* any write. Raises OverrideError."""
    if not isinstance(overrides, list):
        raise OverrideError(
            "overrides payload must be a JSON list of {driver_path, value} objects"
        )
    control_dict_keys, scoped_entries, scopes, regeneration_scopes = _catalog_entries(driver_context)
    scope_by_token = {scope.token: scope for scope in scopes}
    regen_scope_by_key: dict[str, RegenerationScope] = {
        key: regen_scope
        for regen_scope in regeneration_scopes
        for key in regen_scope.selector_keys
    }
    for ov in overrides:
        if not isinstance(ov, dict) or "driver_path" not in ov or "value" not in ov:
            raise OverrideError(
                f"each override must be an object with 'driver_path' and 'value' (got {ov!r})"
            )
        dp = ov["driver_path"]
        if ":" in dp:
            file_path, _, entry_path = dp.partition(":")
            if not _is_safe_system_path(file_path):
                raise OverrideError(f"override file path {file_path!r} is not a safe system/ path")
            if not entry_path:
                raise OverrideError(f"override driver_path {dp!r} is missing an entry path after ':'")
            continue
        elif not dp.startswith("$"):
            if dp in regen_scope_by_key:
                # A bare selector that RESTRUCTURES the file (renames a
                # sub-block, changes which sibling keys are legal) rather
                # than patching one leaf in place -- e.g. myocardiumSolver.
                # Still enum-checked against the catalog exactly like any
                # other entry; only the *application* differs (regenerate
                # vs. key-patch), not the trust boundary.
                entry = scoped_entries.get(dp)
                enum_values = getattr(entry, "enum_values", None)
                if enum_values and ov["value"] not in enum_values:
                    raise OverrideError(
                        f"override {dp!r} value {ov['value']!r} not in enum {tuple(enum_values)}"
                    )
                continue
            # Backward compatibility: flat strings are treated as controlDict entries.
            # Still must be a real controlDict key -- otherwise this silently passes
            # validation and, at apply time, either raises a raw KeyError (no
            # foamDictionary) or silently writes a brand-new bogus key into
            # controlDict (foamDictionary auto-creates missing keys on `-set`).
            if dp not in control_dict_keys:
                known = ", ".join(sorted(control_dict_keys))
                raise OverrideError(
                    f"override driver_path {dp!r} is not a known controlDict entry. "
                    f"Known controlDict entries: {known}"
                )
            continue

        if "<" in dp or ">" in dp:
            raise OverrideError(
                f"override driver_path {dp!r} contains a placeholder; substitute the "
                f"concrete name"
            )

        token = _scope_token(dp)
        if token not in scope_by_token:
            known = ", ".join(f"${t}" for t in sorted(scope_by_token)) or "(none declared)"
            raise OverrideError(
                f"override driver_path {dp!r} uses unknown scope token {'$' + token!r}. "
                f"Known scope tokens: {known}"
            )

        entry = scoped_entries.get(dp)
        if entry is None:
            entry = _match_dynamic_entry(dp, scoped_entries.values())
            if entry is None:
                raise OverrideError(
                    f"override driver_path {dp!r} is not catalog-addressable / applyable"
                )
        enum_values = getattr(entry, "enum_values", None)
        if enum_values and ov["value"] not in enum_values:
            raise OverrideError(
                f"override {dp!r} value {ov['value']!r} not in enum {tuple(enum_values)}"
            )


def apply_overrides(
    overrides: list[dict[str, Any]], *, case_root: Path, driver_context=None,
) -> None:
    """Apply validated overrides to the case dicts.

    Raises OverrideError on any mutator failure (caught at the CLI boundary). Not
    transactional: a mid-list failure can leave earlier overrides applied.
    """
    from openfoam_driver.core.compatibility import resolve_public_driver_context

    driver_context = resolve_public_driver_context(driver_context)
    scope_by_token = {
        scope.token: scope
        for scope in driver_context.capabilities.override_scopes.scopes()
    }
    regen_scope_by_key: dict[str, RegenerationScope] = {
        key: regen_scope
        for regen_scope in driver_context.capabilities.dict_regeneration.scopes()
        for key in regen_scope.selector_keys
    }
    for ov in overrides:
        dp, value = ov["driver_path"], ov["value"]
        try:
            if ":" in dp:
                file_path, _, entry_path = dp.partition(":")
                # foamDictionary spells this scope as "solvers/V/tolerance";
                # update_foam_entry takes it apart. Going through it rather
                # than straight to foamDictionary keeps this route usable
                # without a sourced OpenFOAM, like every other override path.
                *scope_path, key = entry_path.split("/")
                update_foam_entry(
                    case_root / file_path, key, value, scope=scope_path or None
                )
            elif not dp.startswith("$") and dp in regen_scope_by_key:
                regen_scope = regen_scope_by_key[dp]
                # Other $TOKEN. overrides in this same call that target the
                # scope's file: the rebuild needs to see these (not just
                # the selector) because update_foam_entry can only patch a
                # key that already exists, and the new solver may require
                # a key the old file never had (see RegenerationScope.regenerate).
                extra_overrides = {
                    other["driver_path"]: other["value"]
                    for other in overrides
                    if other is not ov
                    and str(other.get("driver_path", "")).startswith("$")
                    and scope_by_token.get(_scope_token(other["driver_path"])) is not None
                    and scope_by_token[_scope_token(other["driver_path"])].file_relpath
                    == regen_scope.file_relpath
                }
                regen_scope.regenerate(
                    case_root / regen_scope.file_relpath, dp, value, extra_overrides,
                )
            elif not dp.startswith("$"):
                update_foam_entry(case_root / "system" / "controlDict", dp, value)
            else:
                token = _scope_token(dp)
                scope = scope_by_token.get(token)
                if scope is None:
                    raise OverrideError(f"unknown scope token {'$' + token!r}")
                scope_path, key = scope.resolve_entry(dp, case_root)
                update_foam_entry(case_root / scope.file_relpath, key, value, scope=scope_path)
        except (OSError, KeyError, ValueError, RuntimeError) as exc:
            raise OverrideError(f"failed to apply override {dp!r}: {exc}") from exc
