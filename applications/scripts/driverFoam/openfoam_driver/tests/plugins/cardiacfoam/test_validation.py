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
#     test_validation
#
# Description
#     Tests validation logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Tests for ``validate_run``.

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
    CONTROL_DICT_ENTRIES,
    DictEntry,
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import ValidationError, slot_key, validate_run

_PHASE_ORDER = ("anatomy", "physics", "stimulus", "solver")


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    yield from CONTROL_DICT_ENTRIES
    for group in get_electro_property_entry_groups().values():
        yield from group


def _blank_run(**overrides) -> RunDocument:
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    for ph, slice_ in overrides.get("config", {}).items():
        config.setdefault(ph, {}).update(slice_)
    return RunDocument(id="r1", name="r", status="draft", config=config)


class TestSlotKeyScopeTokenStripping:
    """slot_key's prefix strip is a syntactic transform (recognize and strip
    a "$SCOPE_TOKEN." shape) -- it must not hardcode the one token the
    built-in cardiac plugin happens to declare, since a future plugin can
    register its own scope token under the same $TOKEN. convention."""

    def test_strips_the_cardiac_scope_token(self) -> None:
        assert slot_key("$ELECTRO_MODEL_COEFFS.myocardiumSolver") == "myocardiumSolver"

    def test_strips_an_arbitrary_scope_token_of_the_same_shape(self) -> None:
        assert slot_key("$SOME_OTHER_PLUGIN_COEFFS.foo.bar") == "foo.bar"

    def test_leaves_an_unprefixed_path_unchanged(self) -> None:
        assert slot_key("myocardiumSolver") == "myocardiumSolver"

    def test_leaves_a_dollar_sign_not_matching_the_scope_token_shape_unchanged(self) -> None:
        # No trailing "." after an all-caps run means this isn't the
        # $TOKEN. convention -- must not be stripped.
        assert slot_key("$notAToken") == "$notAToken"


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
        is_unconditionally_required = e.required and not e.required_when
        is_conditionally_required = e.required_when and any(
            (lambda vals: config.get(ph2, {}).get(k) in (vals if isinstance(vals, tuple) else (vals,)))(v)
            for k, v in e.required_when.items()
            for ph2 in _PHASE_ORDER
        )
        if not (is_unconditionally_required or is_conditionally_required):
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
    # Normalize the stubbed tissue to one compatible with the final ionicModel
    # so the helper yields a genuinely valid run (the generic enum_values[0]
    # fill can otherwise pair, e.g., AlievPanfilov with epicardialCells, which
    # the tissue-compatibility rule now rejects). Unknown models are left as-is.
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    phys = config.get("physics", {})
    model = phys.get("ionicModel")
    if model and "tissue" in phys:
        entry = IONIC_MODEL_CATALOG.get(model)
        if entry and entry.compatible_tissues and phys["tissue"] not in entry.compatible_tissues:
            phys["tissue"] = entry.compatible_tissues[0]
    return RunDocument(id="r1", name="r", status="draft", config=config)


def test_empty_run_reports_missing_required_fields_per_phase():
    errors = validate_run(_blank_run())
    phases_with_errors = {e.phase for e in errors}
    assert {"physics"} <= phases_with_errors
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


# -------- Structured constraint evaluation --------
#
# These tests use synthesized DictEntry fixtures rather than the live
# catalog so the assertions stay stable independently of catalog changes.


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


def test_applicable_when_matches_a_scope_prefixed_predicate_key():
    """applicable_when keys are written in catalog form -- they carry the
    leading "$ELECTRO_MODEL_COEFFS." scope token the same as any other
    driver_path. The predicate lookup must strip that token before
    comparing against context, which is always in slot_key (prefix-
    stripped) form. Regression: before the fix, EVERY applicable_when
    referencing a real driver_path -- not just the bare virtual
    "$..._present"/"$..._supported" tokens -- silently never matched,
    which is what made the restitutionEikonalSolver1D build drop its own
    solver-specific keys without error."""
    from openfoam_driver.specs.dict_builder import select_applicable_entries

    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.gatedByPrefixedKey",
        applicable_when={
            "$ELECTRO_MODEL_COEFFS.verificationModel.type": (
                "manufacturedFDABidomainVerifier",
            ),
        },
    )
    inactive = select_applicable_entries(
        {"verificationModel.type": "manufacturedEikonalVerifier"}, entries=[entry],
    )
    assert inactive == []

    active = select_applicable_entries(
        {"verificationModel.type": "manufacturedFDABidomainVerifier"}, entries=[entry],
    )
    assert active == [entry]


def test_applicable_when_matches_a_dynamic_placeholder_sibling_key():
    """A dynamic_path entry's applicable_when may name a SIBLING leaf inside
    the same <name>-templated block -- e.g. restitutionEikonalSolver1D's
    useEdgeConductance gated on the sibling conductionSystemSolver leaf one
    level up in the same conductionNetworkDomains.<name>.
    purkinjeGraphModelCoeffs block. The predicate key still carries the
    literal "<name>" placeholder while context carries the concrete
    resolved instance name (e.g. "purkinjeNetwork"), so the two can never
    be made equal by prefix-stripping alone -- the placeholder must be
    treated as a wildcard and matched against any configured instance."""
    from openfoam_driver.specs.dict_builder import select_applicable_entries

    entry = _entry(
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>."
        "purkinjeGraphModelCoeffs.useEdgeConductance",
        dynamic_path=True,
        applicable_when={
            "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>."
            "purkinjeGraphModelCoeffs.conductionSystemSolver": (
                "restitutionEikonalSolver1D",
            ),
        },
    )
    inactive = select_applicable_entries(
        {
            "conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs"
            ".conductionSystemSolver": "eikonalSolver1D",
        },
        entries=[entry],
    )
    assert inactive == []

    active = select_applicable_entries(
        {
            "conductionNetworkDomains.purkinjeNetwork.purkinjeGraphModelCoeffs"
            ".conductionSystemSolver": "restitutionEikonalSolver1D",
        },
        entries=[entry],
    )
    assert active == [entry]


def test_validate_run_accepts_default_entries_for_backward_compat():
    """When no entries kwarg is supplied, validate_run uses the live
    catalog (existing public API contract)."""
    run = _filled_run()
    errors = [e for e in validate_run(run) if e.level == "error"]
    assert errors == []


# -------- Solver-coupling evaluator --------
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


def test_solver_coupling_allows_bidomain_with_monodomain1D():
    """bidomainSolver supports monodomain1DSolver via reactionDiffusionPvjCoupler."""
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
    assert len(bidomain_errors) == 0


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


# -------- Block-reference evaluator --------


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


def test_dynamic_required_field_flags_missing_value_scoped_to_its_own_network():
    """purkinjeCV is required_when conductionSystemSolver=eikonalSolver1D,
    but only within the SAME conductionNetworkDomains.<name> block. Two
    networks must be validated independently: a network missing purkinjeCV
    must be flagged even though a sibling network satisfies every
    requirement, and a network that doesn't select eikonalSolver1D must
    never be told it needs purkinjeCV just because another network does."""
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    config["physics"]["myocardiumSolver"] = "eikonalSolver"
    config["physics"][
        "conductionNetworkDomains.networkA.purkinjeGraphModelCoeffs"
        ".conductionSystemSolver"
    ] = "eikonalSolver1D"
    # networkA intentionally omits purkinjeCV.
    config["physics"][
        "conductionNetworkDomains.networkB.purkinjeGraphModelCoeffs"
        ".conductionSystemSolver"
    ] = "monodomain1DSolver"
    # networkB never needs purkinjeCV under this solver.
    run = RunDocument(id="r1", name="r", status="draft", config=config)

    errors = validate_run(run, entries=[])
    cv_errors = [e for e in errors if "purkinjeCV" in e.message]

    assert any("networkA" in e.message for e in cv_errors), (
        f"expected a missing-purkinjeCV error scoped to networkA, got: "
        f"{[e.message for e in errors]}"
    )
    assert not any("networkB" in e.message for e in cv_errors), (
        f"networkB must never be flagged for purkinjeCV -- it does not use "
        f"eikonalSolver1D: {[e.message for e in errors]}"
    )


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
#     test_validation_fixtures
#
# Description
#     Tests validation fixtures logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Cross-fixture regression guard for validate_run.

For each of the 7 tutorial spec fixtures, build a representative RunDocument
that reflects the spec's solver type and assert that ``validate_run`` returns
zero *error*-level violations.

This catches accidentally over-restrictive structured constraints.
Warnings are permitted; only ``level="error"`` must be empty for each fixture run.

Fixture-to-solver mapping (derived from each spec's defaults.ELECTRO_PROPERTIES_SCOPE):
    single_cell          → singleCellSolver
    manufactured_monodomain_pseudo_ecg     → monodomainSolver  (ionic = monodomainFDAManufactured)
    manufactured_bidomain → bidomainSolver (ionic = bidomainFDAManufactured)
    manufactured_bath_bidomain → bidomainSolver (ionic = bathBidomainFDAManufactured)
    niederer_2012        → monodomainSolver  (ionic = TNNP, tissue = epicardialCells)
    restitution_curves   → singleCellSolver  (ionic = TNNP, tissue = epicardialCells)
    generic_case         → monodomainSolver  (representative; generic_case is
                           solver-agnostic, monodomainSolver is the most common)
"""


import pytest

from openfoam_driver.dict_entries import (
    CONTROL_DICT_ENTRIES,
    get_electro_property_entry_groups,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import slot_key, validate_run

_PHASE_ORDER = ("anatomy", "physics", "stimulus", "solver")


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    yield from CONTROL_DICT_ENTRIES
    for group in get_electro_property_entry_groups().values():
        yield from group


def _filled_run_for_solver(myocardium_solver: str, **extra_config) -> RunDocument:
    """Build a RunDocument with every required entry pre-populated.

    Uses the same pattern as ``_filled_run`` in test_validation.py, then
    applies solver-specific overrides so the correct solver is selected and
    any solver-specific required entries are populated.  ``extra_config``
    maps phase → {slot_key: value} for additional overrides.
    """
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    # Pre-populate all required=True entries with plausible stubs.
    # Also pre-populate required_when entries where the predicate matches the
    # known solver so the validator's required_when check passes.
    solver_context = {"myocardiumSolver": myocardium_solver}
    for e in _all_entries():
        is_unconditionally_required = e.required and not e.required_when
        is_conditionally_required = e.required_when and any(
            (lambda vals: solver_context.get(k) in (vals if isinstance(vals, tuple) else (vals,)))(v)
            for k, v in e.required_when.items()
        )
        if not (is_unconditionally_required or is_conditionally_required):
            continue
        ph = next((p for p in _PHASE_ORDER if p in e.phases), None)
        if ph is None:
            continue
        key = slot_key(e.driver_path)
        if e.value_kind == "enum" and e.enum_values:
            config[ph][key] = e.enum_values[0]
        else:
            config[ph][key] = "stub"

    # Apply solver-specific overrides that match the fixture's actual configuration.
    config["physics"]["myocardiumSolver"] = myocardium_solver

    for ph, slice_ in extra_config.items():
        config.setdefault(ph, {}).update(slice_)

    return RunDocument(id="r1", name="r", status="draft", config=config)


# ---------------------------------------------------------------------------
# Fixtures parameterised by spec name + representative run
# ---------------------------------------------------------------------------

# Each tuple is (spec_label, RunDocument).
# The RunDocument is built to match the spec's actual solver and ionic model.

_FIXTURE_RUNS = [
    (
        "single_cell",
        _filled_run_for_solver(
            "singleCellSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "manufactured_monodomain_pseudo_ecg",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "monodomainFDAManufactured",
                # tissue not required for manufactured models (applicable_when excludes them)
            },
        ),
    ),
    (
        "manufactured_bidomain",
        _filled_run_for_solver(
            "bidomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "bidomainFDAManufactured",
            },
        ),
    ),
    (
        "manufactured_bath_bidomain",
        _filled_run_for_solver(
            "bidomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "bathBidomainFDAManufactured",
            },
        ),
    ),
    (
        "niederer_2012",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "restitution_curves",
        _filled_run_for_solver(
            "singleCellSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "generic_case",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
]


@pytest.mark.parametrize("spec_label,run", _FIXTURE_RUNS, ids=[t[0] for t in _FIXTURE_RUNS])
def test_representative_run_has_no_validator_errors(spec_label: str, run: RunDocument):
    """validate_run must return zero error-level violations for each fixture.

    Warnings are permitted.  An error-level violation indicates an over-
    restrictive structured constraint.
    """
    errors = [e for e in validate_run(run) if e.level == "error"]
    assert errors == [], (
        f"spec='{spec_label}': expected no validator errors for representative run, "
        f"got:\n" + "\n".join(f"  [{e.phase}] {e.field}: {e.message}" for e in errors)
    )


# -------- reactionDiffusionPvjCoupler's graph-aware rPvj requirement --------
#
# reactionDiffusionPvjCoupler.C (src/electroModels/electroCouplers/pvjCoupler/
# reactionDiffusion/reactionDiffusionPvjCoupler.C:120-134): if the graph file
# provides per-terminal resistances (a non-empty top-level "pvjResistances"
# list), those are used and a dict-level "rPvj" is never read. Only when the
# graph provides no such list does the C++ side fall back to
# dict.get<scalar>("rPvj") -- a hard FatalError if that key is also absent.
# This is a launch-time semantic check (needs the materialized graph file on
# disk), not a generic catalog rule, so it lives in
# _evaluate_pvj_resistance_requirement and is consumed by the cardiacfoam
# plugin's validate_configuration (the strict pre-flight check gating
# `foamctl run --strict`), not validate_run_semantics.

def _build_pvj_case(tmp_path, *, coupler="reactionDiffusionPvjCoupler",
                     myocardium_solver="monodomainSolver",
                     conduction_solver="monodomain1DSolver",
                     set_rpvj=False, graph_present=None, graph_has_resistances=False,
                     graph_file_key=True):
    from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties

    prefix = (
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork"
        ".purkinjeGraphModelCoeffs"
    )
    overrides = {
        "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.purkinjeNetwork"
        ".conductionSystemDomain": "purkinjeGraphModel",
        f"{prefix}.conductionSystemSolver": conduction_solver,
        f"{prefix}.vm1DRest": "-0.084",
        f"{prefix}.rootStimulus.node": "0",
        f"{prefix}.rootStimulus.startTime": "0.0",
        f"{prefix}.rootStimulus.duration": "0.0",
        f"{prefix}.rootStimulus.intensity": "0.0",
        f"{prefix}.outputVariables.export": "(activationTime)",
        "$ELECTRO_MODEL_COEFFS.domainCouplings.pvj.conductionNetworkDomain": "purkinjeNetwork",
        "$ELECTRO_MODEL_COEFFS.domainCouplings.pvj.couplingMode": "unidirectional",
        "$ELECTRO_MODEL_COEFFS.domainCouplings.pvj.electroDomainCoupler": coupler,
    }
    if conduction_solver == "eikonalSolver1D":
        overrides[f"{prefix}.purkinjeCV"] = "[0 1 -1 0 0 0 0] 4.2"
    if myocardium_solver == "eikonalSolver":
        overrides["$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach"] = "true"
        overrides["$ELECTRO_MODEL_COEFFS.stimulusLocationMin"] = "(1e6 1e6 1e6)"
        overrides["$ELECTRO_MODEL_COEFFS.stimulusLocationMax"] = "(1e6 1e6 1e6)"
    if graph_file_key:
        overrides[f"{prefix}.graphFile"] = "purkinjeGraph"
    if set_rpvj:
        overrides[f"{prefix}.rPvj"] = "150.0"

    selectors = {"myocardiumSolver": myocardium_solver}
    if myocardium_solver != "eikonalSolver":
        selectors["ionicModel"] = "BuenoOrovio"
        selectors["tissue"] = "epicardialCells"
    text = build_electro_properties(selectors, overrides=overrides)
    (tmp_path / "constant").mkdir()
    electro_path = tmp_path / "constant" / "electroProperties"
    electro_path.write_text(text)

    if graph_present:
        graph_text = (
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class dictionary;\n    object purkinjeGraph;\n}\n\n"
            "conductionEdges\n(\n    (0 1 1.0 2.0)\n);\n"
            "pvjNodes (1);\npoints ((0 0 0) (1 0 0));\n"
            "pvjLocations ((1 0 0));\n"
        )
        if graph_has_resistances:
            graph_text += "pvjResistances (150.0);\n"
        (tmp_path / "constant" / "purkinjeGraph").write_text(graph_text)

    return electro_path


def test_pvj_resistance_silent_when_rpvj_explicitly_set(tmp_path):
    """rPvj supplied directly -- valid regardless of graph/resistance state."""
    from openfoam_driver.plugins.cardiacfoam.validation import (
        _evaluate_pvj_resistance_requirement,
    )
    electro_path = _build_pvj_case(tmp_path, set_rpvj=True, graph_present=False)
    diagnostics = _evaluate_pvj_resistance_requirement(tmp_path, electro_path)
    assert diagnostics == ()


def test_pvj_resistance_defers_when_graph_not_yet_materialized(tmp_path):
    """No rPvj, and the referenced graph file does not exist on disk yet
    (cardiacCore generates it later) -- must defer, not error."""
    from openfoam_driver.plugins.cardiacfoam.validation import (
        _evaluate_pvj_resistance_requirement,
    )
    electro_path = _build_pvj_case(tmp_path, set_rpvj=False, graph_present=False)
    diagnostics = _evaluate_pvj_resistance_requirement(tmp_path, electro_path)
    assert diagnostics == ()


def test_pvj_resistance_silent_when_graph_provides_terminal_resistances(tmp_path):
    """No rPvj, but the materialized graph provides a non-empty
    pvjResistances list -- valid, matching reactionDiffusionPvjCoupler.C's
    terminalResistances() precedence over the dict-level rPvj lookup."""
    from openfoam_driver.plugins.cardiacfoam.validation import (
        _evaluate_pvj_resistance_requirement,
    )
    electro_path = _build_pvj_case(
        tmp_path, set_rpvj=False, graph_present=True, graph_has_resistances=True,
    )
    diagnostics = _evaluate_pvj_resistance_requirement(tmp_path, electro_path)
    assert diagnostics == ()


def test_pvj_resistance_errors_when_graph_materialized_without_resistances_and_no_rpvj(tmp_path):
    """No rPvj, graph IS materialized, but it has no pvjResistances -- this
    is the case reactionDiffusionPvjCoupler.C's dict.get<scalar>("rPvj")
    would hard-FatalError on. Neither source exists: must be an error."""
    from openfoam_driver.plugins.cardiacfoam.validation import (
        _evaluate_pvj_resistance_requirement,
    )
    electro_path = _build_pvj_case(
        tmp_path, set_rpvj=False, graph_present=True, graph_has_resistances=False,
    )
    diagnostics = _evaluate_pvj_resistance_requirement(tmp_path, electro_path)
    assert len(diagnostics) == 1
    assert diagnostics[0].level == "error"
    assert "rPvj" in diagnostics[0].message
    assert "purkinjeNetwork" in diagnostics[0].message


def test_pvj_resistance_irrelevant_for_a_different_coupler(tmp_path):
    """The graph-aware rPvj requirement is specific to
    reactionDiffusionPvjCoupler's own dict.get<scalar>("rPvj") fallback --
    eikonalPvjCoupler doesn't read rPvj at all, so this check must never
    fire for it regardless of graph/resistance state."""
    from openfoam_driver.plugins.cardiacfoam.validation import (
        _evaluate_pvj_resistance_requirement,
    )
    electro_path = _build_pvj_case(
        tmp_path, coupler="eikonalPvjCoupler",
        myocardium_solver="eikonalSolver", conduction_solver="eikonalSolver1D",
        set_rpvj=False, graph_present=False,
    )
    diagnostics = _evaluate_pvj_resistance_requirement(tmp_path, electro_path)
    assert diagnostics == ()
