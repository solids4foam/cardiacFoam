import importlib

import pytest

_MODULES = [
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.restitution_curves",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.manufactured_eikonal_ecg",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.manufactured_bath_bidomain",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.manufactured_monodomain_pseudo_ecg",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.cable_1d_cv_convergence",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.cable_1d_restitution",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.manufactured_purkinje_graph",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.manufactured_monodomain_total_lagrangian_em",
    "openfoam_driver.plugins.cardiacfoam.tutorials.defaults.single_cell",
]
_DEAD_NAMES = (
    "POSTPROCESS_SCRIPT_RELPATH",
    "POSTPROCESS_FUNCTION_NAME",
    "CV_EXTRACT_SCRIPT_RELPATH",
    "TABLE_SUMMARY_RELPATH",
)


@pytest.mark.parametrize("module_name", _MODULES)
def test_defaults_module_does_not_export_dead_postprocess_constants(module_name):
    module = importlib.import_module(module_name)
    for name in _DEAD_NAMES:
        assert not hasattr(module, name), f"{module_name} still defines {name}"
        assert name not in getattr(module, "__all__", ()), (
            f"{module_name}.__all__ still lists {name}"
        )
