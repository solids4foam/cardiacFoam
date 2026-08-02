from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_fda_bidomain import (
    make_spec as make_manufactured_fda_bidomain_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_fda_bath_bidomain import (
    make_spec as make_manufactured_fda_bath_bidomain_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_eikonal_ecg import (
    make_spec as make_manufactured_eikonal_ecg_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_monodomain_total_lagrangian_em import (
    make_spec as make_manufactured_monodomain_total_lagrangian_em_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_purkinje_graph import (
    make_spec as make_manufactured_purkinje_graph_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.heart_solver_comparison import (
    make_spec as make_heart_solver_comparison_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.generic_case import make_spec as make_generic_case_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.monodomain_and_eikonal_1d_cable_cv_convergence import (
    make_spec as make_monodomain_and_eikonal_1d_cable_cv_convergence_spec,
)
from openfoam_driver.plugins.cardiacfoam.tutorials.manufactured_fda import make_spec as make_manufactured_fda_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.niederer_2012 import make_spec as make_niederer_2012_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.restitution_curves import make_spec as make_restitution_curves_spec
from openfoam_driver.plugins.cardiacfoam.tutorials.single_cell import make_spec as make_single_cell_spec

SPEC_FACTORIES = {
    "singleCell": make_single_cell_spec,
    "singlecell": make_single_cell_spec,
    "monodomainAndEikonal1DCableCVConvergence": make_monodomain_and_eikonal_1d_cable_cv_convergence_spec,
    "monodomainandeikonal1dcablecvconvergence": make_monodomain_and_eikonal_1d_cable_cv_convergence_spec,
    "niederer2012": make_niederer_2012_spec,
    "niedereretal2012": make_niederer_2012_spec,
    "manufacturedFDA": make_manufactured_fda_spec,
    "manufacturedfda": make_manufactured_fda_spec,
    "manufacturedFDABidomain": make_manufactured_fda_bidomain_spec,
    "manufacturedfdabidomain": make_manufactured_fda_bidomain_spec,
    "manufacturedFDABathBidomain": make_manufactured_fda_bath_bidomain_spec,
    "manufacturedfdabathbidomain": make_manufactured_fda_bath_bidomain_spec,
    "manufacturedEikonalECG": make_manufactured_eikonal_ecg_spec,
    "manufacturedeikonalecg": make_manufactured_eikonal_ecg_spec,
    "manufacturedMonodomainTotalLagrangianEM": make_manufactured_monodomain_total_lagrangian_em_spec,
    "manufacturedelectromechanicsbc": make_manufactured_monodomain_total_lagrangian_em_spec,
    "manufacturedPurkinjeGraph": make_manufactured_purkinje_graph_spec,
    "manufacturedpurkinjegraph": make_manufactured_purkinje_graph_spec,
    "heartSolverComparison": make_heart_solver_comparison_spec,
    "heartsolvercomparison": make_heart_solver_comparison_spec,
    "restitutionCurves": make_restitution_curves_spec,
    "restitutioncurves": make_restitution_curves_spec,
}

REGISTERED_TUTORIALS = (
    "singleCell",
    "monodomainAndEikonal1DCableCVConvergence",
    "niederer2012",
    "manufacturedFDA",
    "manufacturedFDABidomain",
    "manufacturedFDABathBidomain",
    "manufacturedEikonalECG",
    "manufacturedMonodomainTotalLagrangianEM",
    "manufacturedPurkinjeGraph",
    "heartSolverComparison",
    "restitutionCurves",
)
