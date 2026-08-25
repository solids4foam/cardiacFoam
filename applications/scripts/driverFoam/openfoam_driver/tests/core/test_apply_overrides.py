import pytest
from pathlib import Path
from openfoam_driver.tests.conftest import assert_foam_entry, skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.core.specs.apply_overrides import (
    validate_overrides,
    apply_overrides,
    OverrideError,
)
from openfoam_driver.core.runtime.mutators import read_foam_entry

REPO_ROOT = Path(__file__).resolve().parents[6]
SINGLE_CELL = REPO_ROOT / "tutorials" / "electrophysiologyProtocols" / "singleCell"


def _case(tmp_path: Path) -> Path:
    """A minimal case: a real electroProperties + a controlDict carrying deltaT."""
    (tmp_path / "constant").mkdir()
    (tmp_path / "system").mkdir()
    (tmp_path / "system" / "electro").mkdir()
    (tmp_path / "constant" / "electroProperties").write_text(
        (SINGLE_CELL / "constant" / "electroProperties").read_text()
    )
    (tmp_path / "system" / "controlDict").write_text("deltaT    0.001;\nendTime    1;\n")
    (tmp_path / "system" / "fvSolution").write_text("solvers { V { tolerance 1e-5; } }\n")
    (tmp_path / "system" / "electro" / "fvSolution").write_text("solvers { V { tolerance 1e-5; } }\n")
    return tmp_path


# --- validation (catalog-only; no file read) -------------------------------

def test_validate_rejects_non_catalog_driver_path():
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "$ELECTRO_MODEL_COEFFS.notAKey", "value": "1"}])
    assert "notAKey" in str(exc.value)

def test_validate_accepts_safe_file_paths():
    validate_overrides([
        {"driver_path": "system/fvSolution:solvers/V/tolerance", "value": "1e-6"},
        {"driver_path": "system/electro/fvSolution:solvers/V/tolerance", "value": "1e-6"},
        {"driver_path": "system/controlDict:deltaT", "value": "0.0001"},
    ])

def test_validate_rejects_unsafe_file_paths():
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "constant/physicsProperties:type", "value": "electroModel"}])
    assert "not a safe system/ path" in str(exc.value)

    with pytest.raises(OverrideError):
        validate_overrides([{"driver_path": "system/../../etc/passwd:root", "value": "hack"}])

    with pytest.raises(OverrideError):
        validate_overrides([{"driver_path": "/etc/passwd:root", "value": "hack"}])

    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "system/fvSolution:", "value": "1"}])
    assert "missing an entry path" in str(exc.value)


def test_validate_accepts_control_and_electro():
    # flat controlDict, flat electro, nested electro — none should raise.
    validate_overrides([
        {"driver_path": "deltaT", "value": "0.0005"},
        # initialODEStep was removed from the catalog: it is dead (zero reads in
        # this repo or in OpenFOAM). maxSteps is the live ODESolver analogue,
        # read flat from the coeffs dict at ODESolver.C:70.
        {"driver_path": "$ELECTRO_MODEL_COEFFS.maxSteps", "value": "20000"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude", "value": "80"},
    ])


def test_validate_accepts_keys_at_the_scopes_that_actually_read_them():
    """A key catalogued at the top coeffs level is not thereby addressable at
    the nested scopes that separately read it.

    conductionSystemDomain.C:508 binds coeffsDict_ to the
    purkinjeGraphModelCoeffs sub-dict and hands it to ionicModel::New
    (:282), which passes it to ODESolver::New (ionicModel.H:275) -- so
    upstream OpenFOAM reads solver/maxSteps from *that* dict, not from the
    top-level one. Both are used by real tutorials but were rejected by
    validate_overrides, making them unsettable through the CLI.
    """
    validate_overrides([
        {"driver_path": "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.LV"
                        ".purkinjeGraphModelCoeffs.solver", "value": "RKF45"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.LV"
                        ".purkinjeGraphModelCoeffs.maxSteps", "value": "1000"},
    ])


def test_validate_accepts_purkinje_conduction_velocity_and_ode_tolerances():
    """The Purkinje graph carries its own purkinjeCV and its own ODE tolerances.

    purkinjeGraphModelCoeffs.purkinjeCV is a literal conduction velocity in
    m/s (eikonalSolver1D.C:134 -- t = Tact + edgeLength/purkinjeCV), NOT the
    same quantity as the top-level eikonal c0, which has dimensions s^-1/2
    and only becomes a velocity via c0*sqrt(M) (eikonalMyocardiumDomain.C:359).
    The two used to share the bare key 'c0' at different scopes of one file;
    the Purkinje side was renamed to purkinjeCV to remove the ambiguity.

    absTol/relTol reach upstream ODESolver.C:68-69 through the same
    sub-dict binding as the already-catalogued solver/maxSteps.
    """
    purkinje = "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.LV.purkinjeGraphModelCoeffs"
    validate_overrides([
        {"driver_path": f"{purkinje}.purkinjeCV", "value": "[0 1 -1 0 0 0 0] 4.2"},
        {"driver_path": f"{purkinje}.absTol", "value": "1e-8"},
        {"driver_path": f"{purkinje}.relTol", "value": "1e-6"},
    ])


def test_validate_rejects_unknown_flat_controldict_key():
    # A flat (non-$, non-":") driver_path is routed to controlDict for backward
    # compatibility, but must still be a real controlDict key -- not silently
    # accepted and only discovered to be bogus (or worse, silently written) at
    # apply time.
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "notAKey", "value": "1"}])
    assert "notAKey" in str(exc.value)


def test_validate_rejects_malformed_payload():
    with pytest.raises(OverrideError):
        validate_overrides({"driver_path": "deltaT", "value": "1"})   # not a list
    with pytest.raises(OverrideError):
        validate_overrides([{"driver_path": "deltaT"}])               # missing value


def test_validate_rejects_placeholder_path():
    with pytest.raises(OverrideError) as exc:
        validate_overrides(
            [{"driver_path": "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.<AC_name>",
              "value": "1"}]
        )
    assert "contains a placeholder" in str(exc.value)

def test_validate_accepts_concrete_dynamic_path():
    validate_overrides(
        [{"driver_path": "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.myChannel",
          "value": "1"}]
    )


def test_validate_rejects_out_of_enum_value():
    with pytest.raises(OverrideError):
        validate_overrides(
            [{"driver_path": "$ELECTRO_MODEL_COEFFS.solutionAlgorithm", "value": "bogus"}]
        )


# --- application (exercises the real resolver + mutators) ------------------

def test_apply_deltat_edits_control_dict(tmp_path):
    case = _case(tmp_path)
    apply_overrides([{"driver_path": "deltaT", "value": "0.0005"}], case_root=case)
    assert "0.0005" in (case / "system" / "controlDict").read_text()


def test_apply_flat_electro_key_edits_solver_coeffs(tmp_path):
    # maxSteps, not initialODEStep: the latter had no reader anywhere and was
    # deleted from the tutorial dicts, so there is no longer a key to edit.
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "$ELECTRO_MODEL_COEFFS.maxSteps", "value": "20000"}],
        case_root=case,
    )
    assert_foam_entry(
        case / "constant" / "electroProperties",
        "maxSteps",
        "20000",
        scope="singleCellSolverCoeffs",
    )


def test_apply_nested_electro_key_edits_nested_block(tmp_path):
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude", "value": "80"}],
        case_root=case,
    )
    assert "80" in (case / "constant" / "electroProperties").read_text()


def test_apply_missing_control_dict_raises_override_error(tmp_path):
    (tmp_path / "system").mkdir()
    (tmp_path / "constant").mkdir()
    (tmp_path / "constant" / "electroProperties").write_text(
        (SINGLE_CELL / "constant" / "electroProperties").read_text()
    )
    # no controlDict -> mutator FileNotFoundError must surface as OverrideError, not raw.
    with pytest.raises(OverrideError):
        apply_overrides([{"driver_path": "deltaT", "value": "0.0005"}], case_root=tmp_path)


def test_apply_fvSolution_edits_file(tmp_path, monkeypatch):
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/fvSolution:solvers/V/tolerance", "value": "1e-6"}],
        case_root=case,
    )
    fv_solution = case / "system" / "fvSolution"
    tolerance = read_foam_entry(fv_solution, "tolerance", scope=["solvers", "V"])
    assert float(tolerance) == pytest.approx(1e-6)


def test_apply_region_fvSolution_edits_file(tmp_path, monkeypatch):
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/electro/fvSolution:solvers/V/tolerance", "value": "1e-6"}],
        case_root=case,
    )
    electro_fv_solution = case / "system" / "electro" / "fvSolution"
    tolerance = read_foam_entry(electro_fv_solution, "tolerance", scope=["solvers", "V"])
    assert float(tolerance) == pytest.approx(1e-6)
    # Top-level should be unchanged
    top_fv_solution = case / "system" / "fvSolution"
    top_tolerance = read_foam_entry(top_fv_solution, "tolerance", scope=["solvers", "V"])
    assert float(top_tolerance) == pytest.approx(1e-5)


def test_apply_system_file_override_works_without_foamdictionary(tmp_path):
    """A system/<file>:<entry> override must not require a sourced OpenFOAM.

    This route used to call a foamDictionary-only writer directly, so
    fvSolution/fvSchemes overrides only worked in a sourced shell. It now
    goes through the same pure-Python mutator path as the rest of the
    override machinery.
    """
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/fvSolution:solvers/V/tolerance", "value": "1e-9"}],
        case_root=case,
    )
    tolerance = read_foam_entry(
        case / "system" / "fvSolution", "tolerance", scope=["solvers", "V"]
    )
    assert float(tolerance) == pytest.approx(1e-9)


def test_validate_accepts_the_manufactured_solution_switches():
    """The keys that turn manufactured-solution verification on.

    electroVerificationModel::New is called unconditionally
    (myocardiumDomainInterface.C:223); selectedType looks for a
    verificationModel sub-dict and reads `type`, returning nullptr when it is
    absent, empty or "none" (electroVerificationModel.C:42-45). So `type` is
    the one switch -- the redundant `enabled` key it used to share that job
    with is gone, and all four verifier tables now accept "none" alike.

    The Purkinje graph and the PVJ coupling each have their own verifier
    table, so each needs its own `type` at its own scope; the manufactured
    ionic models additionally require `dimension`, which is per-domain (the
    coupled 1D-3D case runs the myocardium at "3D" and the graph at "1D").
    """
    purkinje = "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.LV.purkinjeGraphModelCoeffs"
    validate_overrides([
        {"driver_path": f"{purkinje}.dimension", "value": "1D"},
        {"driver_path": f"{purkinje}.verificationModel.type",
         "value": "manufacturedGraphVerifier"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.domainCouplings.pvj"
                        ".verificationModel.type",
         "value": "coupled1D3DMonodomainVerifier"},
    ])


# --- override scopes are plugin-declared, not core-hardcoded ---------------
#
# P2.5-followup: $ELECTRO_MODEL_COEFFS is no longer the one scope token core
# assumes exists. A plugin declares its own scopes via
# PluginCapabilities.override_scopes; core only knows how to route a "$TOKEN."
# override to whichever scope's token matches.

def test_validate_rejects_an_unknown_scope_token():
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "$SOME_OTHER_PLUGIN.foo", "value": "1"}])
    assert "unknown scope token" in str(exc.value)
    assert "$SOME_OTHER_PLUGIN" in str(exc.value)
    assert "$ELECTRO_MODEL_COEFFS" in str(exc.value)  # the one scope that IS known


def test_apply_rejects_an_unknown_scope_token_before_any_write(tmp_path):
    case = _case(tmp_path)
    with pytest.raises(OverrideError) as exc:
        apply_overrides(
            [{"driver_path": "$SOME_OTHER_PLUGIN.foo", "value": "1"}],
            case_root=case,
        )
    assert "unknown scope token" in str(exc.value)
    # Nothing should have been touched.
    assert "foo" not in (case / "constant" / "electroProperties").read_text()


def test_generic_plugin_declares_zero_override_scopes():
    from openfoam_driver.core.plugin_interface import generic_openfoam_context

    context = generic_openfoam_context()
    assert context.capabilities.override_scopes.scopes() == ()


def test_cardiac_plugin_declares_the_electro_model_coeffs_scope():
    from openfoam_driver.core.plugin_interface import default_driver_context

    context = default_driver_context()
    scopes = context.capabilities.override_scopes.scopes()
    assert len(scopes) == 1
    scope = scopes[0]
    assert scope.token == "ELECTRO_MODEL_COEFFS"
    assert scope.file_relpath == "constant/electroProperties"
    assert scope.catalog_group == "electroProperties"


def test_cardiac_scope_resolve_entry_matches_the_old_hardcoded_behavior(tmp_path):
    """Same (scope_path, key) shape apply_overrides used to compute inline
    via detect_myocardium_solver_name + _entry_scope_and_key."""
    from openfoam_driver.core.plugin_interface import default_driver_context

    case = _case(tmp_path)
    context = default_driver_context()
    scope = context.capabilities.override_scopes.scopes()[0]
    scope_path, key = scope.resolve_entry(
        "$ELECTRO_MODEL_COEFFS.solutionAlgorithm", case,
    )
    assert scope_path == ["singleCellSolverCoeffs"]
    assert key == "solutionAlgorithm"


# --- regeneration scopes: myocardiumSolver routed through override channel -
#
# myocardiumSolver is a bare (non-"$") catalog entry whose value change
# RESTRUCTURES electroProperties (renames the active <solver>Coeffs
# sub-block, flips which sibling keys the catalog allows) rather than
# patching one leaf in place. A plugin declares this via
# PluginCapabilities.dict_regeneration, the sibling of override_scopes for
# $TOKEN. leaves.

PURKINJE_NIEDERER = REPO_ROOT / "tutorials" / "NiedererEtAl2011" / "purkinjeNiedererEtAl2011"


def test_generic_plugin_declares_zero_regeneration_scopes():
    from openfoam_driver.core.plugin_interface import generic_openfoam_context

    context = generic_openfoam_context()
    assert context.capabilities.dict_regeneration.scopes() == ()


def test_cardiac_plugin_declares_the_myocardium_solver_regeneration_scope():
    from openfoam_driver.core.plugin_interface import default_driver_context

    context = default_driver_context()
    scopes = context.capabilities.dict_regeneration.scopes()
    assert len(scopes) == 1
    scope = scopes[0]
    assert scope.selector_keys == frozenset({"myocardiumSolver"})
    assert scope.file_relpath == "constant/electroProperties"
    assert scope.catalog_group == "electroProperties"


def test_validate_accepts_a_bare_myocardium_solver_override_with_a_valid_enum_value():
    # Previously rejected outright: "not a known controlDict entry" (VERIFIED
    # FACT #1 in the task this test guards). Now routed to regeneration.
    validate_overrides([{"driver_path": "myocardiumSolver", "value": "eikonalSolver"}])


def test_validate_rejects_a_bare_myocardium_solver_override_with_an_invalid_enum_value():
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "myocardiumSolver", "value": "notARealSolver"}])
    assert "not in enum" in str(exc.value)


def test_validate_still_rejects_other_bare_selector_keys_as_unknown_controlDict_entries():
    """ionicModel/tissue/conductivitySource are also _SELECTOR_KEYS, but are
    NOT wired into regeneration (see dict_builder._SELECTOR_KEYS docstring
    and electro_properties_regeneration_scope): they change a value in
    place without renaming anything, so they stay on the ordinary
    $ELECTRO_MODEL_COEFFS.* key-patch route. A bare (un-scoped) override
    for one of them is still just an unrecognised controlDict entry."""
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "ionicModel", "value": "TNNP"}])
    assert "not a known controlDict entry" in str(exc.value)


def test_apply_regenerates_electro_properties_for_a_myocardium_solver_override(tmp_path):
    """End-to-end through the public apply_overrides() entry point (not the
    dict_builder function directly): on a copy of the real
    purkinjeNiedererEtAl2011 monodomain fixture, myocardiumSolver=eikonalSolver
    -- bundled, in the same call, with the handful of new fields eikonalSolver
    requires with no catalog default -- must now be REJECTED, propagated as
    an OverrideError.

    The fixture's Purkinje network uses conductionSystemSolver=
    monodomain1DSolver / electroDomainCoupler=reactionDiffusionPvjCoupler.
    regenerate_electro_properties carries that network forward verbatim
    rather than resynthesising it, and eikonalSolver+monodomain1DSolver is
    explicitly invalid per solver_coupling.SOLVER_COMPATIBILITY_RULES
    ("eikonal myocardium cannot couple to reaction-diffusion Purkinje") --
    confirmed against the real, hand-authored electroProperties.eikonal
    fixture for this same tutorial, which uses a different, compatible
    pairing (eikonalSolver1D / eikonalPvjCoupler) instead of a bare carry-
    forward. See test_dict_builder.py's
    test_purkinje_niederer_monodomain_to_eikonal_end_to_end for the same
    scenario exercised directly against regenerate_electro_properties."""
    if not (PURKINJE_NIEDERER / "constant" / "electroProperties.monodomain").exists():
        pytest.skip("tutorial fixture not present in this checkout")

    (tmp_path / "constant").mkdir()
    original_text = (
        PURKINJE_NIEDERER / "constant" / "electroProperties.monodomain"
    ).read_text()
    (tmp_path / "constant" / "electroProperties").write_text(original_text)

    overrides = [
        {"driver_path": "myocardiumSolver", "value": "eikonalSolver"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach", "value": "true"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.stimulusLocationMin", "value": "(1e6 1e6 1e6)"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.stimulusLocationMax", "value": "(1e6 1e6 1e6)"},
    ]
    validate_overrides(overrides)
    with pytest.raises(OverrideError) as exc:
        apply_overrides(overrides, case_root=tmp_path)
    message = str(exc.value)
    assert "eikonal" in message.lower()
    assert "reaction-diffusion" in message.lower()

    # The original file must be untouched -- regenerate_electro_properties
    # stages the rewrite and only commits it once the final, carried-forward
    # result validates.
    assert (tmp_path / "constant" / "electroProperties").read_text() == original_text
