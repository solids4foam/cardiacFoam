"""Focused guards for repository contracts repaired after the consistency audit."""

from pathlib import Path

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.plugins.cardiacfoam.active_tension_catalog import ACTIVE_TENSION_MODEL_CATALOG


REPO_ROOT = Path(__file__).resolve().parents[7]


def test_manufactured_active_tension_is_introspectable():
    entry = ACTIVE_TENSION_MODEL_CATALOG["ManufacturedElectromechanics"]
    assert entry.states == ("Ta",)
    assert entry.algebraic == (
        "AV_Vm",
        "AV_lambda",
        "AV_voltageActivation",
        "AV_lengthFactor",
    )
    assert entry.constants == ("AC_Tmax", "AC_V0", "AC_gamma")


def test_dashboard_is_not_advertised_by_package_contract():
    pyproject = (REPO_ROOT / "applications/scripts/driverFoam/pyproject.toml").read_text()
    readme = (
        REPO_ROOT / "applications/scripts/driverFoam/openfoam_driver/README.md"
    ).read_text()
    assert "dashboard = [" not in pyproject
    assert "driverFoam dashboard" not in readme
    assert "openfoam_driver.dashboard" not in readme


def test_readme_uses_registered_compact_batched_selectors():
    readme = (REPO_ROOT / "src/ionicModels/README.md").read_text()
    for model in (
        "AlievPanfilov",
        "BuenoOrovio",
        "Courtemanche",
        "Fabbri",
        "Gaur",
        "Grandi",
        "PerisYague",
        "Stewart",
        "TNNP",
        "ToRORd_dynCl",
        "Trovato",
        "TWorld",
    ):
        assert f"`{model}compactBatched`" in readme


def test_eikonal_ecg_rejects_non_transmural_modes_before_weighting():
    source = (
        REPO_ROOT / "src/electroModels/ecgModels/eikonalECG/eikonalECG.C"
    ).read_text()
    guard = 'if (mode != "transmuralBands")'
    weighting = "ionicHeterogeneity::transmuralBandWeights"
    assert guard in source
    assert source.index(guard) < source.index(weighting)
    assert "namedRegions" not in source
    assert "cellZoneRegions" not in source
