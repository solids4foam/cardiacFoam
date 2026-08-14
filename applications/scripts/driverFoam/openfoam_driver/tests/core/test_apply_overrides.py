import os

import pytest
from pathlib import Path
from openfoam_driver.tests.conftest import assert_foam_entry, skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.specs.apply_overrides import (
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
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "$ELECTRO_MODEL_COEFFS.initialODEStep", "value": "2e-5"}],
        case_root=case,
    )
    assert_foam_entry(
        case / "constant" / "electroProperties",
        "initialODEStep",
        "2e-5",
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
    import shutil
    if not shutil.which("foamDictionary"):
        if os.environ.get("REQUIRE_FOAMDICTIONARY") == "1":
            pytest.fail("foamDictionary is required for this integration test")
        pytest.skip("foamDictionary not available")
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/fvSolution:solvers/V/tolerance", "value": "1e-6"}],
        case_root=case,
    )
    fv_solution = case / "system" / "fvSolution"
    tolerance = read_foam_entry(fv_solution, "tolerance", scope=["solvers", "V"])
    assert float(tolerance) == pytest.approx(1e-6)


def test_apply_region_fvSolution_edits_file(tmp_path, monkeypatch):
    import shutil
    if not shutil.which("foamDictionary"):
        if os.environ.get("REQUIRE_FOAMDICTIONARY") == "1":
            pytest.fail("foamDictionary is required for this integration test")
        pytest.skip("foamDictionary not available")
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
