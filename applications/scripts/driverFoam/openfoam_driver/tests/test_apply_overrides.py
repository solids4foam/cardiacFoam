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
    (tmp_path / "constant" / "electroProperties").write_text(
        (SINGLE_CELL / "constant" / "electroProperties").read_text()
    )
    (tmp_path / "system" / "controlDict").write_text("deltaT    0.001;\nendTime    1;\n")
    return tmp_path


# --- validation (catalog-only; no file read) -------------------------------

def test_validate_rejects_non_catalog_driver_path():
    with pytest.raises(OverrideError) as exc:
        validate_overrides([{"driver_path": "notAKey", "value": "1"}])
    assert "notAKey" in str(exc.value)


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
    with pytest.raises(OverrideError):
        validate_overrides(
            [{"driver_path": "$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.<AC_name>",
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
