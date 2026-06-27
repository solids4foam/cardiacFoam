import pytest
from pathlib import Path

from openfoam_driver.specs.apply_overrides import (
    validate_overrides,
    apply_overrides,
    OverrideError,
)

REPO_ROOT = Path(__file__).resolve().parents[5]
SINGLE_CELL = REPO_ROOT / "tutorials" / "singleCellprotocols" / "singleCell"


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
        {"driver_path": "$ELECTRO_MODEL_COEFFS.initialODEStep", "value": "2e-5"},
        {"driver_path": "$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude", "value": "80"},
    ])


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
    text = (case / "constant" / "electroProperties").read_text()
    assert "initialODEStep" in text and "2e-5" in text


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
        pytest.skip("foamDictionary not available")
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/fvSolution:solvers/V/tolerance", "value": "1e-6"}],
        case_root=case,
    )
    assert "1e-6" in (case / "system" / "fvSolution").read_text()


def test_apply_region_fvSolution_edits_file(tmp_path, monkeypatch):
    import shutil
    if not shutil.which("foamDictionary"):
        pytest.skip("foamDictionary not available")
    case = _case(tmp_path)
    apply_overrides(
        [{"driver_path": "system/electro/fvSolution:solvers/V/tolerance", "value": "1e-6"}],
        case_root=case,
    )
    assert "1e-6" in (case / "system" / "electro" / "fvSolution").read_text()
    assert "1e-5" in (case / "system" / "fvSolution").read_text()  # Top-level should be unchanged
