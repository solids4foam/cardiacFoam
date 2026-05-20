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


# -------- P5e.1/2: solver-coupling evaluator --------
#
# The three prose-only entries (conductionSystemSolver, electroDomainCoupler,
# conductionNetworkDomain) don't fit the four DictEntry families, but the
# rules they encode are already machine-readable via
# SOLVER_COMPATIBILITY_RULES in solver_coupling.py. These tests pin the
# behaviour we expect from _evaluate_solver_coupling.


def _coupling_run(myocardium: str, *,
                  purkinje: str | None = None,
                  coupler: str | None = None,
                  network_name: str = "purkinjeNet",
                  coupling_name: str = "lvCoupling") -> RunDocument:
    """Build a run with selected solver + optional Purkinje pairing.

    Dynamic-path slot_keys (e.g. domainCouplings.lvCoupling.electroDomainCoupler)
    are written into the physics slice — matches how _flatten_context will
    expose them.
    """
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    config["physics"]["myocardiumSolver"] = myocardium
    if purkinje is not None:
        config["physics"][
            f"conductionNetworkDomains.{network_name}."
            f"purkinjeGraphModelCoeffs.conductionSystemSolver"
        ] = purkinje
        # The network must be declared as a block — i.e. at least one
        # sub-key exists under conductionNetworkDomains.<name>.*.
        config["physics"][
            f"conductionNetworkDomains.{network_name}.purkinjeGraphModelCoeffs.someKey"
        ] = "x"
    if coupler is not None:
        config["physics"][
            f"domainCouplings.{coupling_name}.electroDomainCoupler"
        ] = coupler
        # A coupling references the network by name.
        config["physics"][
            f"domainCouplings.{coupling_name}.conductionNetworkDomain"
        ] = network_name
    return RunDocument(id="r1", name="r", status="draft", config=config)


def test_solver_coupling_silent_when_no_purkinje_pairing():
    """No conductionSystemSolver in context → no coupling rules fire."""
    run = _coupling_run("monodomainSolver")
    errors = validate_run(run, entries=[])
    coupling_errors = [
        e for e in errors
        if "coupling" in e.message.lower() or "coupler" in e.message.lower()
    ]
    assert coupling_errors == []


def test_solver_coupling_valid_monodomain_pair_silent():
    """Valid pair (mono + monodomain1D + reactionDiffusionPvjCoupler)
    must not emit any solver-coupling error."""
    run = _coupling_run(
        "monodomainSolver",
        purkinje="monodomain1DSolver",
        coupler="reactionDiffusionPvjCoupler",
    )
    errors = validate_run(run, entries=[])
    coupling_errors = [
        e for e in errors
        if "incompatible" in e.message.lower()
        or "required_coupler" in e.message.lower()
    ]
    assert coupling_errors == [], (
        f"valid pair must not error, got: {[e.message for e in errors]}"
    )


def test_solver_coupling_flags_incompatible_mono_eikonal_pair():
    """mono myocardium + eikonal Purkinje is invalid per the rules table."""
    run = _coupling_run(
        "monodomainSolver",
        purkinje="eikonalSolver",
        coupler="reactionDiffusionPvjCoupler",
    )
    errors = validate_run(run, entries=[])
    incompat = [e for e in errors if "incompatible" in e.message.lower()]
    assert len(incompat) >= 1, (
        f"expected incompatible-pair error, got: {[e.message for e in errors]}"
    )


def test_solver_coupling_flags_bidomain_with_purkinje():
    """bidomainSolver does not support Purkinje coupling at all."""
    run = _coupling_run(
        "bidomainSolver",
        purkinje="monodomain1DSolver",
        coupler="reactionDiffusionPvjCoupler",
    )
    errors = validate_run(run, entries=[])
    bidomain_errors = [
        e for e in errors
        if "bidomain" in e.message.lower() and "purkinje" in e.message.lower()
    ]
    assert len(bidomain_errors) >= 1


def test_solver_coupling_flags_wrong_coupler_for_valid_pair():
    """Valid mono+monodomain1D pair but the wrong coupler → error citing
    the required_coupler."""
    run = _coupling_run(
        "monodomainSolver",
        purkinje="monodomain1DSolver",
        coupler="eikonalPvjCoupler",   # wrong; should be reactionDiffusionPvjCoupler
    )
    errors = validate_run(run, entries=[])
    coupler_errors = [
        e for e in errors
        if "reactiondiffusionpvjcoupler" in e.message.lower()
    ]
    assert len(coupler_errors) >= 1, (
        f"expected error citing reactionDiffusionPvjCoupler, got: "
        f"{[e.message for e in errors]}"
    )


# -------- P5e.3/4: block-reference evaluator --------


def test_block_reference_silent_when_no_couplings():
    """No domainCouplings in context → no block-reference rules fire."""
    run = _coupling_run("monodomainSolver")
    errors = validate_run(run, entries=[])
    ref_errors = [e for e in errors if "reference" in e.message.lower()]
    assert ref_errors == []


def test_block_reference_silent_when_target_block_declared():
    """conductionNetworkDomain references a name that has at least one
    sub-key under conductionNetworkDomains.<name>.* → no error."""
    run = _coupling_run(
        "monodomainSolver",
        purkinje="monodomain1DSolver",
        coupler="reactionDiffusionPvjCoupler",
        network_name="purkinjeNet",
        coupling_name="lvCoupling",
    )
    errors = validate_run(run, entries=[])
    dangling_errors = [
        e for e in errors
        if "reference" in e.message.lower()
        and ("not declared" in e.message.lower() or "dangling" in e.message.lower())
    ]
    assert dangling_errors == []


def test_block_reference_flags_dangling_target():
    """conductionNetworkDomain points at a name that has no matching block
    declaration → error."""
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    config["physics"]["myocardiumSolver"] = "monodomainSolver"
    config["physics"][
        "domainCouplings.lvCoupling.conductionNetworkDomain"
    ] = "ghostNet"   # never declared under conductionNetworkDomains.ghostNet.*
    run = RunDocument(id="r1", name="r", status="draft", config=config)

    errors = validate_run(run, entries=[])
    dangling = [
        e for e in errors
        if "ghostNet" in e.message
        and ("not declared" in e.message.lower()
             or "no matching" in e.message.lower())
    ]
    assert len(dangling) >= 1, (
        f"expected dangling-reference error for ghostNet, got: "
        f"{[e.message for e in errors]}"
    )
