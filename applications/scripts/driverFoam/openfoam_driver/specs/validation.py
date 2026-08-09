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
#     validation
#
# Description
#     Evaluates simulation constraints and flags missing required fields.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Run-document validator.

The validator reports three kinds of issue:

1. **Required-field omissions.** Each ``DictEntry`` flagged ``required=True``
   must have a value in the Run document slice owned by its *primary* phase.
2. **Enum violations.** Entries with ``value_kind="enum"`` whose value is not
   one of the declared ``enum_values``.
3. **Structured constraints.** Each ``DictEntry`` may declare
   ``applicable_when``, ``forbidden_when``, ``required_when``, and
   ``mutually_exclusive_with``. The validator evaluates these against a
   flattened view of the run config. The legacy hardcoded ``eikonalSolver``/
   ``ionicModel`` cross-field check is now encoded as
   ``forbidden_when={"myocardiumSolver": "eikonalSolver"}`` on the
   ``ionicModel`` entry and evaluated programmatically here.

The *primary* phase is the first phase in workflow order
(``anatomy → physics → stimulus → solver``) that the entry claims.
Multi-phase entries are validated there; the other phases do not duplicate
validation errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Iterable

from openfoam_driver.core.contracts.dictionary import DictEntry, Phase

if TYPE_CHECKING:
    from openfoam_driver.core.plugin_interface import DriverContext

_PHASE_ORDER: tuple[Phase, ...] = (
    "anatomy", "physics", "stimulus", "solver",
)


from .validation_types import ValidationError
def _all_entries(driver_context: "DriverContext | None" = None):
    if driver_context is None:
        from openfoam_driver.core.plugin_interface import default_driver_context
        driver_context = default_driver_context()
    yield from driver_context.plugin.get_dict_entries()


def primary_phase(entry) -> Phase | None:
    """Return the editing phase for a (possibly multi-phase) entry.

    Walks ``_PHASE_ORDER`` and returns the first phase the entry claims;
    every other declared phase is a read-only mirror.
    """
    for ph in _PHASE_ORDER:
        if ph in entry.phases:
            return ph
    return None


_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."


def slot_key(driver_path: str) -> str:
    """Map a driver_path to its slot key inside a phase slice.

    Strips the ``$ELECTRO_MODEL_COEFFS.`` prefix when present; otherwise
    returns the path as-is. Multi-segment unprefixed paths are kept intact
    so that nested-group leaves don't collide with top-level keys of the
    same name (e.g. ``$ELECTRO_MODEL_COEFFS.bathPotentialDomain.phiEReferenceValue`` must not overwrite the
    top-level ``type`` entry inside the physics slice).
    """
    if driver_path.startswith(_COEFFS_PREFIX):
        return driver_path[len(_COEFFS_PREFIX):]
    return driver_path


def _slice_value(run, phase: Phase, driver_path: str):
    """Look up the slot value for a driver_path inside a phase slice."""
    slice_ = run.config.get(phase, {}) or {}
    return slice_.get(slot_key(driver_path))


def _flatten_context(run) -> dict[str, Any]:
    """Build a flat predicate-key → value view across every phase slice.

    Structured-constraint predicates reference slot-keys (the post-
    ``slot_key`` form: ``myocardiumSolver``, ``ionicModel``,
    ``singleCellStimulus.stim_amplitude``, ...). Multiple phases never
    write the same slot-key in practice; if they do, the last wins —
    document the convention rather than silently merging.
    """
    context: dict[str, Any] = {}
    for slice_ in run.config.values():
        if not slice_:
            continue
        for key, val in slice_.items():
            if val not in (None, ""):
                context[key] = val
    return context


def _predicate_matches(
    context: dict[str, Any],
    key: str,
    expected: str | tuple[str, ...],
) -> bool:
    """Return True iff ``context[key]`` matches ``expected``.

    Scalar ``expected`` → equality. Tuple ``expected`` → membership.
    A missing key is treated as not-matching (the predicate's
    precondition is absent).
    """
    if key not in context:
        return False
    actual = _normalise_word(context[key])
    if isinstance(expected, tuple):
        return actual in tuple(_normalise_word(item) for item in expected)
    return actual == _normalise_word(expected)


def _normalise_word(value: Any) -> Any:
    """Strip one balanced OpenFOAM word/string quote pair for comparisons."""
    if not isinstance(value, str) or len(value) < 2:
        return value
    if (value[0], value[-1]) in {('"', '"'), ("'", "'")}:
        return value[1:-1]
    return value


def _entry_is_applicable(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Evaluate ``applicable_when`` and ``forbidden_when`` constraints.

    Returns False if any ``forbidden_when`` predicate matches. Otherwise,
    returns True if all ``applicable_when`` predicates match (or if
    ``applicable_when`` is empty).
    """
    if entry.forbidden_when:
        if any(
            _predicate_matches(context, key, expected)
            for key, expected in entry.forbidden_when.items()
        ):
            return False

    if not entry.applicable_when:
        return True
    return all(
        _predicate_matches(context, key, expected)
        for key, expected in entry.applicable_when.items()
    )


def is_required_in_context(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Whether an entry's ``required`` semantics fire under this context.

    Many entries declare BOTH ``required=True`` and ``required_when={...}``
    — the author's intent is "required, but only when the predicate
    matches". This helper reads the two fields together:

    - ``required_when`` non-empty → required iff any predicate matches.
    - ``required_when`` empty     → ``entry.required`` is taken at face value.
    """
    if entry.required_when:
        return any(
            _predicate_matches(context, key, expected)
            for key, expected in entry.required_when.items()
        )
    return entry.required


def _entry_value_present(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Is the entry's own slot set in the flattened context?"""
    key = slot_key(entry.driver_path)
    return key in context and context[key] not in (None, "")


def _format_predicate(predicate: dict[str, Any]) -> str:
    parts: list[str] = []
    for key, expected in predicate.items():
        if isinstance(expected, tuple):
            parts.append(f"{key} ∈ {{{', '.join(expected)}}}")
        else:
            parts.append(f"{key}={expected}")
    return " and ".join(parts)


def _all_entries_list(driver_context: "DriverContext | None" = None):
    return list(_all_entries(driver_context))


def validate_run(
    run,
    *,
    entries: Iterable[DictEntry] | None = None,
    driver_context: "DriverContext | None" = None,
) -> list[ValidationError]:
    """Validate ``run`` against the dict-entry catalog.

    ``entries`` overrides the live catalog for testability and for callers
    that want to validate against a curated subset (e.g., dict_builder).
    When omitted, the full live catalog is used.
    """
    if driver_context is None:
        from openfoam_driver.core.plugin_interface import default_driver_context
        driver_context = default_driver_context()
    entry_list: list[DictEntry] = (
        list(entries) if entries is not None else _all_entries_list(driver_context)
    )
    context = _flatten_context(run)
    errors: list[ValidationError] = []

    # 1) Required-field checks. Skip entries whose applicable_when fails —
    #    requiredness is conditional on applicability. When an entry has
    #    both ``required=True`` and a non-empty ``required_when``, the
    #    latter narrows the former: requiredness fires only when at least
    #    one ``required_when`` predicate matches.
    for e in entry_list:
        if not _entry_is_applicable(e, context):
            continue
        if not is_required_in_context(e, context):
            continue
        if e.dynamic_path:
            # Dynamic-path entries describe templates; concrete required
            # leaves are the user's responsibility when those blocks are
            # actually configured.
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val in (None, ""):
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{e.driver_path} is required.",
                level="error",
            ))

    # 2) Enum checks. Skip inapplicable entries for the same reason.
    for e in entry_list:
        if e.value_kind != "enum" or not e.enum_values:
            continue
        if not _entry_is_applicable(e, context):
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val is None or val == "":
            continue
        normalised_val = _normalise_word(val)
        normalised_enum_values = tuple(_normalise_word(item) for item in e.enum_values)
        if normalised_val not in normalised_enum_values:
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{val!r} is not one of {list(e.enum_values)}.",
                level="error",
            ))

    # 3) Structured constraints.
    # (The ionicModel entry carries forbidden_when={"myocardiumSolver": "eikonalSolver"}
    # which the section below evaluates programmatically.)
    errors.extend(_evaluate_structured(entry_list, context))

    # 4) Domain semantics are a plugin concern.  Core owns only generic
    # catalog constraints and receives solver-specific diagnostics as data.
    errors.extend(driver_context.plugin.validate_run_semantics(context))

    return errors


def _evaluate_structured(
    entries: list[DictEntry],
    context: dict[str, Any],
) -> list[ValidationError]:
    """Evaluate the four structured-constraint families per entry."""
    errors: list[ValidationError] = []
    paths_set = {slot_key(e.driver_path) for e in entries
                 if _entry_value_present(e, context)}

    for e in entries:
        ph = primary_phase(e) or "physics"

        # forbidden_when: fires when ANY predicate matches AND the entry's
        # own slot has a value. Each matching predicate emits its own
        # ValidationError so the reason text stays specific.
        if _entry_value_present(e, context):
            for key, expected in e.forbidden_when.items():
                if _predicate_matches(context, key, expected):
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is forbidden when "
                            f"{_format_predicate({key: expected})}."
                        ),
                        level="error",
                    ))

        # Skip entries whose applicable_when/forbidden_when precondition fails
        if not _entry_is_applicable(e, context):
            continue

        # required_when: handled by section 1 of validate_run via
        # is_required_in_context. Section 3 does NOT re-emit a violation
        # to avoid double-firing on the same entry. The structured
        # required_when field is still consumed — its predicates feed
        # is_required_in_context which gates section 1's required check.

        # mutually_exclusive_with: fires when BOTH this entry's slot is set
        # AND any of the listed sibling slots is also set. To avoid
        # double-reporting symmetric relations we only flag the side that
        # *declares* the relation.
        if _entry_value_present(e, context):
            for sibling_path in e.mutually_exclusive_with:
                sibling_slot = slot_key(sibling_path)
                if sibling_slot in paths_set:
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is mutually exclusive with "
                            f"{sibling_path}."
                        ),
                        level="error",
                    ))

    return errors

