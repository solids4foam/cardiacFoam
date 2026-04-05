from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class DictEntry:
    driver_path: str
    description: str
    source_refs: tuple[str, ...]
    notes: str = ""
    value_kind: str = "openfoam_literal"
    ui_control: str = "text"
    enum_values: tuple[str, ...] = ()
    examples: tuple[str, ...] = ()
    dynamic_path: bool = False


PHYSICS_PROPERTY_ENTRIES: Final[tuple[DictEntry, ...]] = (
    DictEntry(
        driver_path="type",
        description=(
            "Top-level physics model selector. Cardiac tutorial values in this repository "
            "include electroModel and electroMechanicalModel."
        ),
        source_refs=(
            "modules/physicsModel/src/solids4FoamModels/physicsModel/physicsModel.C",
            "applications/utilities/listCellModelsVariables/listCellModelsVariables.C",
            "src/electroModels/core/electroModel/electroModel.H",
        ),
        value_kind="enum",
        ui_control="select",
        enum_values=("electroModel", "electroMechanicalModel"),
    ),
)


ELECTRO_PROPERTY_ENTRY_GROUPS: Final[dict[str, tuple[DictEntry, ...]]] = {
    "top_level": (
        DictEntry(
            driver_path="myocardiumSolver",
            description="Myocardium diffusion solver. Selects the numerical kernel for 3D reaction-diffusion. Reads from 'myocardiumSolverCoeffs' sub-dict.",
            source_refs=(
                "src/electroModels/domains/myocardiumDomain/reactionDiffusionSolver.C",
                "src/electroModels/core/electroModel/electrophysicsSystemBuilder.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("monodomainSolver", "bidomainSolver", "eikonalSolver"),
        ),
    ),
    "common_model_coeffs": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
            description="Electro solver time-discretisation mode.",
            source_refs=(
                "src/electroModels/core/electroModel/electroModel.C",
                "src/verificationModels/manufacturedFDA/monodomainVerification/manufacturedFDAMonodomainVerifier.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("implicit", "explicit"),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicModel",
            description="Ionic cell model selector.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/electroModels/wrappers/electroMechanicalModel/electroMechanicalModel.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=(
                "AlievPanfilov",
                "BuenoOrovio",
                "Courtemanche",
                "Fabbri",
                "Gaur",
                "Grandi",
                "ORd",
                "Stewart",
                "TNNP",
                "ToRORd_dynCl",
                "Trovato",
                "tmanufacturedFDA",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.tissue",
            description="Tissue selector for ionic-model specialisation.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicSelector.C",
                "src/electroModels/wrappers/electroMechanicalModel/electroMechanicalModel.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("epicardialCells", "mCells", "endocardialCells", "myocyte"),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.dimension",
            description="Dimensional selector used by manufactured/verification ionic models.",
            source_refs=("src/ionicModels/ionicModel/ionicSelector.C",),
            value_kind="enum",
            ui_control="select",
            enum_values=("1D", "2D", "3D"),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.writeAfterTime",
            description="Suppresses single-cell trace output before the given time.",
            source_refs=("src/genericWriter/ionicModelIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.utilities",
            description="Utility mode flag for ionic-model helper applications.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="boolean",
            ui_control="checkbox",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initSampleCell",
            description="Integration-point index used for sampled single-cell style output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="integer",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.export",
            description="Ionic variables exported to volumetric fields or trace output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="word_list",
            ui_control="token_list",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.debug",
            description="Ionic variables printed in debug output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="word_list",
            ui_control="token_list",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.export",
            description="Active-tension variables exported to fields or traces.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            ui_control="token_list",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.debug",
            description="Active-tension variables printed in debug output.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            ui_control="token_list",
        ),
    ),
    "ode_solver_passthrough": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solver",
            description="ODE solver selector passed through to OpenFOAM's ODESolver factory.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.H",
                "src/activeTensionModels/GoktepeKuhl/GoktepeKuhl.C",
                "src/activeTensionModels/NashPanfilov/NashPanfilov.C",
            ),
            notes=(
                "The repository source shows pass-through to ODESolver::New(*this, dict_). "
                "Additional ODESolver-specific keys may exist beyond the commonly used entries "
                "listed here."
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("RKF45",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initialODEStep",
            description="Initial ODE step-size hint passed through to the selected ODE solver.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.maxSteps",
            description="Maximum internal ODE steps allowed per macro time step.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="integer",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.absTol",
            description="Absolute tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.relTol",
            description="Relative tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
        ),
    ),
    "single_cell_stimulus": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start",
            description="Start time for the S1/S2 stimulus protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1",
            description="S1 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_duration",
            description="Pulse duration for both S1 and S2 stimuli.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
            description="Stimulus amplitude for the S1/S2 protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim1",
            description="Number of S1 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S2",
            description="S2 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim2",
            description="Number of S2 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            ui_control="number",
        ),
    ),
    "monodomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivity",
            description="Monodomain conductivity tensor.",
            source_refs=(
                "src/electroModels/solvers/monoDomainSolver/monoDomainDiffusionSolver.C",
                "src/electroModels/solvers/eikonalSolver/eikonalSolver.C",
            ),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.chi",
            description="Surface-to-volume ratio.",
            source_refs=(
                "src/electroModels/solvers/monoDomainSolver/monoDomainDiffusionSolver.C",
                "src/electroModels/solvers/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.Cm",
            description="Membrane capacitance.",
            source_refs=(
                "src/electroModels/solvers/monoDomainSolver/monoDomainDiffusionSolver.C",
                "src/electroModels/solvers/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.infoFrequency",
            description="Logging cadence for monodomain explicit stepping.",
            source_refs=("src/electroModels/solvers/monoDomainSolver/monoDomainSolver.C",),
            value_kind="integer",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.allowIonicStimulusInMonodomain",
            description="Explicit opt-in for combining ionic-model and PDE stimulus definitions.",
            source_refs=("src/electroModels/solvers/monoDomainSolver/monoDomainSolver.C",),
            value_kind="boolean",
            ui_control="checkbox",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusLocationMin",
            description="Single-box minimum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            ui_control="vector3",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusLocationMax",
            description="Single-box maximum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            ui_control="vector3",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusLocationMinList",
            description="Per-box minimum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusLocationMaxList",
            description="Per-box maximum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusStartTime",
            description="Shared start time for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusStartTimeList",
            description="Per-box start times for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusDuration",
            description="Shared duration for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            ui_control="text",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusDurationList",
            description="Per-box durations for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusIntensity",
            description="Shared current density for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            ui_control="text",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.monodomainStimulus.stimulusIntensityList",
            description="Per-box current densities for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solverHookFields.preProcess",
            description="Manufactured-FDA pre-process field list.",
            source_refs=("src/verificationModels/manufacturedFDA/monodomainVerification/manufacturedFDAMonodomainVerifier.C",),
            value_kind="word_list",
            ui_control="token_list",
        ),
    ),
    "eikonal_diffusion": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach",
            description="Toggle between the advection-diffusion and deferred-correction forms.",
            source_refs=("src/electroModels/solvers/eikonalSolver/eikonalSolver.C",),
            value_kind="boolean",
            ui_control="checkbox",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMin",
            description="Minimum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/solvers/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            ui_control="vector3",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMax",
            description="Maximum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/solvers/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            ui_control="vector3",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.c0",
            description="Wave-speed parameter for the eikonal formulation.",
            source_refs=("src/electroModels/solvers/eikonalDiffusionElectro/eikonalDiffusionElectro.C",),
            value_kind="scalar",
            ui_control="number",
        ),
    ),
    "ecg": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ECGDomain",
            description="ECG domain model selector within the ecgDomains sub-dictionary.",
            source_refs=(
                "src/electroModels/ecgModels/pseudoECGDomain/pseudoECGDomain.H",
                "src/electroModels/core/ecgDomain/ecgDomain.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("pseudoECGDomain",),
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.pseudoECGDomainCoeffs.electrodes.<electrode>",
            description="Per-electrode position vector in the pseudoECGDomain coefficients block.",
            source_refs=("src/electroModels/ecgModels/pseudoECGDomain/pseudoECGDomain.C",),
            notes="The driver can update existing electrode entries but does not insert new electrode names.",
            value_kind="vector3",
            ui_control="vector3",
            dynamic_path=True,
        ),
    ),
    "bidomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityIntracellular",
            description="Intracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/solvers/biDomainSolver/biDomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityExtracellular",
            description="Extracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/solvers/biDomainSolver/biDomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceCell",
            description="Cell index used to pin the extracellular potential reference.",
            source_refs=("src/electroModels/solvers/biDomainSolver/biDomainSolver.C",),
            value_kind="integer",
            ui_control="number",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceValue",
            description="Value of the extracellular potential at the reference cell.",
            source_refs=("src/electroModels/solvers/biDomainSolver/biDomainSolver.C",),
            value_kind="scalar",
            ui_control="number",
        ),
    ),
    "conduction_system": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.conductionSystemDomainModel",
            description=(
                "Conduction system domain model selector (sub-dictionary within "
                "conductionNetworkDomains). Used exclusively in Purkinje/1D pipelines."
            ),
            source_refs=("src/electroModels/conductionSystemModels/purkinjeNetworkModel/purkinjeNetworkModel.H",),
            value_kind="enum",
            ui_control="select",
            enum_values=("purkinjeNetworkModel",),
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCouplingModel",
            description=(
                "PVJ coupling model selector within the domainCouplings sub-dictionary. "
                "Used to couple a Purkinje network domain to the 3-D myocardium."
            ),
            source_refs=("src/electroModels/couplers/pvjResistanceCoupler.H",),
            value_kind="enum",
            ui_control="select",
            enum_values=("pvjResistanceCoupler",),
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.R_pvj",
            description="PVJ coupling resistance [Ohm·m²].",
            source_refs=("src/electroModels/couplers/pvjResistanceCoupler.H",),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjRadius",
            description="Sphere radius for current distribution at PVJ junctions [m].",
            source_refs=("src/electroModels/couplers/pvjResistanceCoupler.H",),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
            description="PVJ coupling direction: networkToMyocardium (1-way) or bidirectional.",
            source_refs=("src/electroModels/couplers/pvjResistanceCoupler.H",),
            value_kind="enum",
            ui_control="select",
            enum_values=("networkToMyocardium", "bidirectional"),
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.advanceScheme",
            description=(
                "Electrophysics time-advance scheme selector. "
                "staggeredElectrophysicsAdvanceScheme staggers the ionic ODE and PDE solves."
            ),
            source_refs=("src/electroModels/core/electroModel/schemes/staggered/staggeredElectrophysicsAdvanceScheme.H",),
            value_kind="enum",
            ui_control="select",
            enum_values=("staggeredElectrophysicsAdvanceScheme",),
        ),
    ),
    "active_tension": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.activeTensionModel",
            description="Active-tension model selector.",
            source_refs=(
                "src/electroModels/wrappers/electroMechanicalModel/electroMechanicalModel.C",
                "src/activeTensionModels/activeTensionModel/activeTensionModel.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("GoktepeKuhl", "NashPanfilov"),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.couplingSignal",
            description="Coupling signal requested by the active-tension model.",
            source_refs=(
                "src/activeTensionModels/GoktepeKuhl/GoktepeKuhl.C",
                "src/activeTensionModels/NashPanfilov/NashPanfilov.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("Vm",),
        ),
    ),
}


def all_documented_driver_paths() -> tuple[str, ...]:
    paths = [entry.driver_path for entry in PHYSICS_PROPERTY_ENTRIES]
    for entries in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        paths.extend(entry.driver_path for entry in entries)
    return tuple(dict.fromkeys(paths))
