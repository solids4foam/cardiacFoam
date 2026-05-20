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
    DictEntry,
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


# -------- P5c: structured constraint evaluation --------
#
# These tests use synthesized DictEntry fixtures rather than the live
# catalog so the assertions stay stable as P5b migrates real entries.


def _entry(driver_path: str, **overrides) -> DictEntry:
    """Tiny DictEntry builder for structured-constraint tests."""
    defaults = {
        "driver_path": driver_path,
        "description": "fixture",
        "source_refs": ("ref.C",),
        "phases": frozenset({"physics"}),
    }
    defaults.update(overrides)
    return DictEntry(**defaults)


def test_forbidden_when_flags_violation_in_run():
    """forbidden_when matches AND entry value is set in context → error."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.ionicModel",
        forbidden_when={"myocardiumSolver": "eikonalSolver"},
    )
    run = _blank_run(config={"physics": {
        "myocardiumSolver": "eikonalSolver",
        "ionicModel": "TNNP",
    }})
    errors = validate_run(run, entries=[entry])
    forbidden_errors = [e for e in errors if "forbidden" in e.message.lower()]
    assert len(forbidden_errors) == 1, (
        f"expected exactly one forbidden_when violation, got: "
        f"{[e.message for e in errors]}"
    )
    assert "myocardiumSolver" in forbidden_errors[0].message
    assert "eikonalSolver" in forbidden_errors[0].message


def test_forbidden_when_silent_when_predicate_doesnt_match():
    """Same entry; non-matching context → no forbidden_when violation."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.ionicModel",
        forbidden_when={"myocardiumSolver": "eikonalSolver"},
    )
    run = _blank_run(config={"physics": {
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
    }})
    errors = validate_run(run, entries=[entry])
    forbidden_errors = [e for e in errors if "forbidden" in e.message.lower()]
    assert forbidden_errors == []


def test_required_when_flags_missing_value():
    """required_when matches AND entry value missing → error."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
        required_when={"myocardiumSolver": "singleCellSolver"},
        phases=frozenset({"stimulus"}),
    )
    run = _blank_run(config={"physics": {
        "myocardiumSolver": "singleCellSolver",
    }})
    errors = validate_run(run, entries=[entry])
    required_errors = [
        e for e in errors
        if "required" in e.message.lower() and "stim_amplitude" in e.message
    ]
    assert len(required_errors) == 1, (
        f"expected one required_when violation, got: {[e.message for e in errors]}"
    )


def test_required_when_silent_when_value_present():
    """required_when matches AND value present → no violation."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
        required_when={"myocardiumSolver": "singleCellSolver"},
        phases=frozenset({"stimulus"}),
    )
    run = _blank_run(config={
        "physics": {"myocardiumSolver": "singleCellSolver"},
        "stimulus": {"singleCellStimulus.stim_amplitude": "60"},
    })
    errors = validate_run(run, entries=[entry])
    assert errors == []


def test_required_when_silent_when_predicate_doesnt_match():
    """required_when doesn't fire when its context predicate doesn't match."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
        required_when={"myocardiumSolver": "singleCellSolver"},
        phases=frozenset({"stimulus"}),
    )
    run = _blank_run(config={"physics": {
        "myocardiumSolver": "monodomainSolver",
    }})
    errors = validate_run(run, entries=[entry])
    assert errors == []


def test_applicable_when_skips_inapplicable_entry():
    """An entry whose applicable_when predicate fails should be entirely
    skipped — even required=True does not fire."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.bidomainOnlyKey",
        required=True,
        applicable_when={"myocardiumSolver": "bidomainSolver"},
    )
    run = _blank_run(config={"physics": {
        "myocardiumSolver": "monodomainSolver",
    }})
    errors = validate_run(run, entries=[entry])
    assert errors == [], (
        f"inapplicable entry must not fire required check, got: "
        f"{[e.message for e in errors]}"
    )


def test_mutually_exclusive_with_flags_violation():
    """When both this entry and a mutex sibling are set → error.

    mutually_exclusive_with paths must be unambiguous — full driver_path
    or slot_key form. Leaf-only names are not supported (collision risk
    across nested groups).
    """
    entry_a = _entry(
        "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDuration",
        mutually_exclusive_with=(
            "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
        ),
        phases=frozenset({"stimulus"}),
    )
    entry_b = _entry(
        "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
        phases=frozenset({"stimulus"}),
    )
    run = _blank_run(config={"stimulus": {
        "externalStimulus.stimulusDuration": "0.002",
        "externalStimulus.stimulusDurationList": "(0.002 0.001)",
    }})
    errors = validate_run(run, entries=[entry_a, entry_b])
    mutex_errors = [e for e in errors if "mutually exclusive" in e.message.lower()]
    assert len(mutex_errors) >= 1, (
        f"expected mutually-exclusive violation, got: {[e.message for e in errors]}"
    )


def test_tuple_predicate_matches_membership():
    """A tuple-valued predicate is satisfied by membership."""
    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.manufacturedCoeff",
        applicable_when={"ionicModel": (
            "monodomainFDAManufactured",
            "bidomainFDAManufactured",
            "bathBidomainFDAManufactured",
        )},
        required=True,
    )
    # Not applicable → no required-check fire.
    run_inactive = _blank_run(config={"physics": {"ionicModel": "TNNP"}})
    assert validate_run(run_inactive, entries=[entry]) == []

    # Applicable → required fires when value missing.
    run_active = _blank_run(config={"physics": {
        "ionicModel": "monodomainFDAManufactured",
    }})
    errors = validate_run(run_active, entries=[entry])
    required_errors = [e for e in errors if "required" in e.message.lower()]
    assert len(required_errors) >= 1


def test_validate_run_accepts_default_entries_for_backward_compat():
    """When no entries kwarg is supplied, validate_run uses the live
    catalog (existing public API contract)."""
    run = _filled_run()
    errors = [e for e in validate_run(run) if e.level == "error"]
    assert errors == []
