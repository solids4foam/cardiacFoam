"""Tests for ``validate_run`` (Phase A Task A6).

The validator walks every ``DictEntry``, attributes errors to the entry's
*primary* (editing) phase, and reports three kinds of issue: required-field
omissions, enum violations, and the small set of v1 cross-field
constraints (currently: ``eikonalSolver`` is incompatible with an explicit
``ionicModel``).

Test 2 uses ``_filled_run`` rather than the plan's hand-written minimal
config because ``dict_entries.py`` declares many required leaf-name keys
inside ``$ELECTRO_MODEL_COEFFS.*``; supplying them programmatically keeps
the test focused on validator behaviour rather than field enumeration.
"""

from __future__ import annotations

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import ValidationError, slot_key, validate_run

_PHASE_ORDER = ("anatomy", "physics", "stimulus", "solver")


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        yield from group


def _blank_run(**overrides) -> RunDocument:
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    for ph, slice_ in overrides.get("config", {}).items():
        config.setdefault(ph, {}).update(slice_)
    return RunDocument(id="r1", name="r", status="draft", config=config)


def _filled_run(**overrides) -> RunDocument:
    """A Run with every required leaf-name pre-populated with a plausible stub.

    Keeps the validator's required-field check happy so the test can isolate
    the behaviour we care about (no errors when the run is complete; or a
    constraint violation when the user toggles an incompatible combination).
    """
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    for e in _all_entries():
        if not e.required:
            continue
        ph = next((p for p in _PHASE_ORDER if p in e.phases), None)
        if ph is None:
            continue
        config.setdefault(ph, {})
        key = slot_key(e.driver_path)
        if e.value_kind == "enum" and e.enum_values:
            config[ph][key] = e.enum_values[0]
        else:
            config[ph][key] = "stub"
    for ph, slice_ in overrides.get("config", {}).items():
        config.setdefault(ph, {}).update(slice_)
    return RunDocument(id="r1", name="r", status="draft", config=config)


def test_empty_run_reports_missing_required_fields_per_phase():
    errors = validate_run(_blank_run())
    phases_with_errors = {e.phase for e in errors}
    assert {"physics", "solver"} <= phases_with_errors
    assert all(isinstance(e, ValidationError) for e in errors)


def test_valid_minimal_run_has_no_errors():
    run = _filled_run()
    errors = [e for e in validate_run(run) if e.level == "error"]
    assert errors == [], f"expected no errors, got: {errors}"


def test_constraint_violation_is_flagged():
    # eikonalSolver disallows an explicit ionicModel.
    run = _filled_run(config={
        "physics": {
            "myocardiumSolver": "eikonalSolver",
            "ionicModel": "tenTusscher2006",
        },
    })
    errors = validate_run(run)
    assert any("eikonal" in e.message.lower() for e in errors), (
        f"expected an eikonal-related error, got: {[e.message for e in errors]}"
    )
