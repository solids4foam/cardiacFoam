from enum import Enum


class CardiacTutorialID(str, Enum):
    SINGLE_CELL = "singleCell"
    CABLE_1D_CV_CONVERGENCE = "cable1DCVConvergence"
    NIEDERER_2012 = "niederer2012"
    MANUFACTURED_MONODOMAIN_PSEUDO_ECG = "manufacturedMonodomainPseudoECG"
    MANUFACTURED_BIDOMAIN = "manufacturedBidomain"
    MANUFACTURED_BATH_BIDOMAIN = "manufacturedBathBidomain"
    MANUFACTURED_EIKONAL_ECG = "manufacturedEikonalECG"
    MANUFACTURED_MONODOMAIN_TOTAL_LAGRANGIAN_EM = "manufacturedMonodomainTotalLagrangianEM"
    MANUFACTURED_MONODOMAIN_1D3D = "manufacturedMonodomain1D3D"
    MANUFACTURED_PURKINJE_GRAPH = "manufacturedPurkinjeGraph"
    HEART_SOLVER_COMPARISON = "heartSolverComparison"
    RESTITUTION_CURVES = "restitutionCurves"
    CABLE_1D_RESTITUTION = "cable1DRestitution"
