from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_bidomain import (
    make_spec as make_manufactured_bidomain_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_bath_bidomain import (
    make_spec as make_manufactured_bath_bidomain_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_eikonal_ecg import (
    make_spec as make_manufactured_eikonal_ecg_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_total_lagrangian_em import (
    make_spec as make_manufactured_monodomain_total_lagrangian_em_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_1d3d import (
    make_spec as make_manufactured_monodomain_1d3d_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_purkinje_graph import (
    make_spec as make_manufactured_purkinje_graph_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.heart_solver_comparison import (
    make_spec as make_heart_solver_comparison_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.generic_case import make_spec as make_generic_case_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.cable_1d_cv_convergence import (
    make_spec as make_cable_1d_cv_convergence_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_pseudo_ecg import make_spec as make_manufactured_monodomain_pseudo_ecg_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.niederer_2012 import make_spec as make_niederer_2012_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.restitution_curves import make_spec as make_restitution_curves_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.single_cell import make_spec as make_single_cell_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.cable_1d_restitution import (
    make_spec as make_cable_1d_restitution_spec,
)

from openfoam_driver.plugins.cardiacfoam.tutorials.ids import CardiacTutorialID

SPEC_FACTORIES = {
    CardiacTutorialID.SINGLE_CELL.value: make_single_cell_spec,
    CardiacTutorialID.SINGLE_CELL.value.lower(): make_single_cell_spec,
    CardiacTutorialID.CABLE_1D_CV_CONVERGENCE.value: make_cable_1d_cv_convergence_spec,
    CardiacTutorialID.CABLE_1D_CV_CONVERGENCE.value.lower(): make_cable_1d_cv_convergence_spec,
    CardiacTutorialID.NIEDERER_2012.value: make_niederer_2012_spec,
    CardiacTutorialID.NIEDERER_2012.value.lower(): make_niederer_2012_spec,
    "niedereretal2012": make_niederer_2012_spec,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_PSEUDO_ECG.value: make_manufactured_monodomain_pseudo_ecg_spec,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_PSEUDO_ECG.value.lower(): make_manufactured_monodomain_pseudo_ecg_spec,
    CardiacTutorialID.MANUFACTURED_BIDOMAIN.value: make_manufactured_bidomain_spec,
    CardiacTutorialID.MANUFACTURED_BIDOMAIN.value.lower(): make_manufactured_bidomain_spec,
    CardiacTutorialID.MANUFACTURED_BATH_BIDOMAIN.value: make_manufactured_bath_bidomain_spec,
    CardiacTutorialID.MANUFACTURED_BATH_BIDOMAIN.value.lower(): make_manufactured_bath_bidomain_spec,
    CardiacTutorialID.MANUFACTURED_EIKONAL_ECG.value: make_manufactured_eikonal_ecg_spec,
    CardiacTutorialID.MANUFACTURED_EIKONAL_ECG.value.lower(): make_manufactured_eikonal_ecg_spec,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_TOTAL_LAGRANGIAN_EM.value: make_manufactured_monodomain_total_lagrangian_em_spec,
    "manufacturedelectromechanicsbc": make_manufactured_monodomain_total_lagrangian_em_spec,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_1D3D.value: make_manufactured_monodomain_1d3d_spec,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_1D3D.value.lower(): make_manufactured_monodomain_1d3d_spec,
    CardiacTutorialID.MANUFACTURED_PURKINJE_GRAPH.value: make_manufactured_purkinje_graph_spec,
    CardiacTutorialID.MANUFACTURED_PURKINJE_GRAPH.value.lower(): make_manufactured_purkinje_graph_spec,
    CardiacTutorialID.HEART_SOLVER_COMPARISON.value: make_heart_solver_comparison_spec,
    CardiacTutorialID.HEART_SOLVER_COMPARISON.value.lower(): make_heart_solver_comparison_spec,
    CardiacTutorialID.RESTITUTION_CURVES.value: make_restitution_curves_spec,
    CardiacTutorialID.RESTITUTION_CURVES.value.lower(): make_restitution_curves_spec,
    CardiacTutorialID.CABLE_1D_RESTITUTION.value: make_cable_1d_restitution_spec,
    CardiacTutorialID.CABLE_1D_RESTITUTION.value.lower(): make_cable_1d_restitution_spec,
}

REGISTERED_TUTORIALS = (
    CardiacTutorialID.SINGLE_CELL.value,
    CardiacTutorialID.CABLE_1D_CV_CONVERGENCE.value,
    CardiacTutorialID.NIEDERER_2012.value,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_PSEUDO_ECG.value,
    CardiacTutorialID.MANUFACTURED_BIDOMAIN.value,
    CardiacTutorialID.MANUFACTURED_BATH_BIDOMAIN.value,
    CardiacTutorialID.MANUFACTURED_EIKONAL_ECG.value,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_TOTAL_LAGRANGIAN_EM.value,
    CardiacTutorialID.MANUFACTURED_MONODOMAIN_1D3D.value,
    CardiacTutorialID.MANUFACTURED_PURKINJE_GRAPH.value,
    CardiacTutorialID.HEART_SOLVER_COMPARISON.value,
    CardiacTutorialID.RESTITUTION_CURVES.value,
    CardiacTutorialID.CABLE_1D_RESTITUTION.value,
)
