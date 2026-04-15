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
    required: bool = False
    constraints: tuple[str, ...] = ()


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
            "src/electroModels/core/electroModel.H",
        ),
        value_kind="enum",
        ui_control="select",
        enum_values=("electroModel", "electroMechanicalModel"),
        required=True,
    ),
)


ELECTRO_PROPERTY_ENTRY_GROUPS: Final[dict[str, tuple[DictEntry, ...]]] = {
    "top_level": (
        DictEntry(
            driver_path="myocardiumSolver",
            description="Top-level myocardium solver selector. Determines the active '<solver>Coeffs' sub-dictionary in electroProperties.",
            source_refs=(
                "src/electroModels/core/electroModel.C",
                "src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("monodomainSolver", "bidomainSolver", "singleCellSolver", "eikonalSolver"),
            required=True,
        ),
    ),
    "common_model_coeffs": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
            description="Electro solver time-discretisation mode.",
            source_refs=(
                "src/electroModels/core/electroModel.C",
                "src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("implicit", "explicit"),
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicModel",
            description="Ionic cell model selector.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
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
                "monodomainFDAManufactured",
                "bidomainFDAManufactured",
            ),
            required=True,
            constraints=("Not applicable when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.tissue",
            description="Tissue selector for ionic-model specialisation.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicSelector.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("epicardialCells", "mCells", "endocardialCells", "myocyte"),
            required=True,
            constraints=("Not applicable when myocardiumSolver=eikonalSolver or ionicModel is a manufactured model.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme",
            description="Time-advance scheme for multi-domain coupling (myocardium, Purkinje, ECG).",
            source_refs=(
                "src/electroModels/core/advanceSchemes/electrophysicsAdvanceScheme.H",
                "src/electroModels/core/advanceSchemes/staggered/staggeredElectrophysicsAdvanceScheme.C",
                "src/electroModels/core/advanceSchemes/pimpleStaggered/pimpleStaggeredElectrophysicsAdvanceScheme.C",
            ),
            notes="staggeredElectrophysicsAdvanceScheme: weakly coupled, fast, stable for unidirectional. pimpleStaggeredElectrophysicsAdvanceScheme: strongly coupled with PIMPLE iteration, stable for bidirectional coupling (requires solutionAlgorithm=implicit).",
            value_kind="enum",
            ui_control="select",
            enum_values=("staggeredElectrophysicsAdvanceScheme", "pimpleStaggeredElectrophysicsAdvanceScheme"),
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.dimension",
            description="Dimensional selector used by manufactured/verification ionic models.",
            source_refs=("src/ionicModels/ionicModel/ionicSelector.C",),
            value_kind="enum",
            ui_control="select",
            enum_values=("1D", "2D", "3D"),
            required=False,
            constraints=("Only applicable for manufactured ionic models (monodomainFDAManufactured, bidomainFDAManufactured).",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.writeAfterTime",
            description="Suppresses single-cell trace output before the given time.",
            source_refs=("src/genericWriter/ionicModelIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.utilities",
            description="Utility mode flag for ionic-model helper applications.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="boolean",
            ui_control="checkbox",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initSampleCell",
            description="Integration-point index used for sampled single-cell style output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="integer",
            ui_control="number",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.export",
            description="Ionic variables exported to volumetric fields or trace output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="word_list",
            ui_control="token_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.debug",
            description="Ionic variables printed in debug output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="word_list",
            ui_control="token_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.export",
            description="Active-tension variables exported to fields or traces.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            ui_control="token_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.debug",
            description="Active-tension variables printed in debug output.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            ui_control="token_list",
            required=False,
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
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initialODEStep",
            description="Initial ODE step-size hint passed through to the selected ODE solver.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.maxSteps",
            description="Maximum internal ODE steps allowed per macro time step.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="integer",
            ui_control="number",
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.absTol",
            description="Absolute tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.relTol",
            description="Relative tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            ui_control="number",
            required=False,
        ),
    ),
    "single_cell_stimulus": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start",
            description="Start time for the S1/S2 stimulus protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1",
            description="S1 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_duration",
            description="Pulse duration for both S1 and S2 stimuli.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
            description="Stimulus amplitude for the S1/S2 protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim1",
            description="Number of S1 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S2",
            description="S2 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim2",
            description="Number of S2 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
    ),
    "monodomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivity",
            description="Monodomain conductivity tensor.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
            required=True,
            constraints=("Required for monodomainSolver and eikonalSolver; not used by singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.chi",
            description="Surface-to-volume ratio.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required for monodomainSolver and bidomainSolver; not used by singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.cm",
            description="Membrane capacitance.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required for monodomainSolver and bidomainSolver; not used by singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.infoFrequency",
            description="Logging cadence for monodomain explicit stepping.",
            source_refs=("src/electroModels/myocardiumModels/monodomainSolver/monoDomainSolver.C",),
            value_kind="integer",
            ui_control="number",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.allowIonicStimulusInMonodomain",
            description="Explicit opt-in for combining ionic-model and PDE stimulus definitions.",
            source_refs=("src/electroModels/myocardiumModels/monodomainSolver/monoDomainSolver.C",),
            value_kind="boolean",
            ui_control="checkbox",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMin",
            description="Single-box minimum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            ui_control="vector3",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMinList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMax",
            description="Single-box maximum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            ui_control="vector3",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMaxList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMinList",
            description="Per-box minimum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            ui_control="textarea",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMin.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMaxList",
            description="Per-box maximum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            ui_control="textarea",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMax.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTime",
            description="Shared start time for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            ui_control="number",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTimeList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTimeList",
            description="Per-box start times for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTime.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDuration",
            description="Shared duration for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            ui_control="text",
            required=False,
            constraints=("Mutually exclusive with stimulusDurationList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
            description="Per-box durations for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
            required=False,
            constraints=("Mutually exclusive with stimulusDuration.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity",
            description="Shared current density for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            ui_control="text",
            required=False,
            constraints=("Mutually exclusive with stimulusIntensityList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensityList",
            description="Per-box current densities for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            ui_control="textarea",
            required=False,
            constraints=("Mutually exclusive with stimulusIntensity.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solverHookFields.preProcess",
            description="Manufactured-FDA pre-process field list.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",),
            value_kind="word_list",
            ui_control="token_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.type",
            description="Optional myocardium-side verification hook selector.",
            source_refs=(
                "src/verificationModels/electroVerification/electroVerificationModel.C",
                "src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.H",
                "src/verificationModels/bidomainVerification/manufacturedFDABidomainVerifier.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=(
                "manufacturedFDAMonodomainVerifier",
                "manufacturedFDABidomainVerifier",
            ),
            required=False,
        ),
    ),
    "eikonal_diffusion": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach",
            description="Toggle between the advection-diffusion and deferred-correction forms.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="boolean",
            ui_control="checkbox",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMin",
            description="Minimum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            ui_control="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMax",
            description="Maximum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            ui_control="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.c0",
            description="Wave-speed parameter for the eikonal formulation.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="scalar",
            ui_control="number",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
    ),
    "ecg": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            description="ECG solver selector within the ecgDomains sub-dictionary.",
            source_refs=(
                "src/electroModels/electroDomains/ecgDomain/ecgSolver.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("pseudoECG",),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.enabled",
            description="Enable manufactured pseudo-ECG verification for the selected ECG domain.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="boolean",
            ui_control="checkbox",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.dimension",
            description="Dimensional selector used by the manufactured pseudo-ECG reference.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="enum",
            ui_control="select",
            enum_values=("1D", "2D", "3D"),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.referenceQuadratureOrder",
            description="Quadrature order used for the manufactured pseudo-ECG reference integral.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="integer",
            ui_control="number",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.checkQuadratureOrders",
            description="Additional quadrature orders used to compare manufactured pseudo-ECG reference convergence.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="label_list",
            ui_control="textarea",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.electrodePositions.<electrode>",
            description="Per-electrode position vector [m] in the ECG domain block.",
            source_refs=("src/electroModels/electroDomains/ecgDomain/ecgDomain.C",),
            notes="The driver can update existing electrode entries but does not insert new electrode names.",
            value_kind="vector3",
            ui_control="vector3",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
    ),
    "bidomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityIntracellular",
            description="Intracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
            required=True,
            constraints=("Required for bidomainSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityExtracellular",
            description="Extracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            ui_control="textarea",
            required=True,
            constraints=("Required for bidomainSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceCell",
            description="Cell index used to pin the extracellular potential reference.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="integer",
            ui_control="number",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceValue",
            description="Value of the extracellular potential at the reference cell.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="scalar",
            ui_control="number",
            required=False,
        ),
    ),
    "conduction_system": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.conductionSystemDomain",
            description=(
                "Conduction system domain model selector (sub-dictionary within "
                "conductionNetworkDomains). Used exclusively in Purkinje/1D pipelines."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("purkinjeNetworkModel",),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when conductionNetworkDomains block is present.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.graphFile",
            description=(
                "REQUIRED: Name of the graph file in constant/ directory "
                "(e.g., 'purkinjeGraph'). File must contain 'edges', 'points', "
                "'pvjNodes', 'pvjLocations' dictionaries/lists."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word",
            ui_control="text",
            dynamic_path=True,
            required=True,
            constraints=("Required within a conductionNetworkDomains entry.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.rootNode",
            description=(
                "Optional: index of the Purkinje root node (default: 0). "
                "Must be a valid node in the graph (0 <= rootNode < nNodes)."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="label",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.conductionSystemSolver",
            description=(
                "1D graph solver used within the conduction network domain. "
                "Default is monodomain1DSolver."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemSolver.C",
                "src/electroModels/conductionSystemModels/monodomain1DSolver/monodomain1DSolver.H",
                "src/electroModels/conductionSystemModels/eikonalSolver1D/eikonalSolver1D.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("monodomain1DSolver", "eikonalSolver"),
            dynamic_path=True,
            required=False,
            constraints=("monodomain1DSolver valid only with monodomainSolver myocardium; eikonalSolver valid only with eikonalSolver myocardium.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.ionicModel",
            description=(
                "Ionic model used for 1D Purkinje monodomain. "
                "Uses the same runtime selection table as the 3D myocardium solver. "
                "Stewart (human Purkinje) is the canonical choice for conduction system simulations."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
                "src/ionicModels/ionicModel/ionicModel.C",
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
                "monodomainFDAManufactured",
                "bidomainFDAManufactured",
            ),
            notes="Stewart is the canonical human Purkinje model. monodomainFDAManufactured/bidomainFDAManufactured are for verification only.",
            dynamic_path=True,
            required=False,
            constraints=("Required when conductionSystemSolver=monodomain1DSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.tissue",
            description=(
                "Tissue type for 1D Purkinje ionic model. "
                "Examples: epicardialCells, endocardialCells, mCells, myocyte."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
                "src/ionicModels/ionicModel/ionicModel.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("epicardialCells", "endocardialCells", "mCells", "myocyte"),
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.pvjNodes",
            description=(
                "Label list of Purkinje graph node indices that are PVJ terminals. "
                "Must match pvjLocations in length."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="label_list",
            ui_control="textarea",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.pvjLocations",
            description=(
                "List of 3D coordinates [m] for each PVJ terminal, one per pvjNodes entry."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="point_list",
            ui_control="textarea",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.pvjRadius",
            description=(
                "Sphere radius [m] around each PVJ centroid used to identify the 3D cells "
                "that receive the coupling current. Default: 0.5e-3."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.rootStimulus.startTime",
            description="Start time [s] of the root-node stimulus applied to Purkinje node 0.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.rootStimulus.duration",
            description="Duration [s] of the root-node stimulus pulse.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.rootStimulus.intensity",
            description="Amplitude [A/m³] of the root-node stimulus current applied to Purkinje node 0.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.chi",
            description="Surface-to-volume ratio [1/m] for the 1D Purkinje monodomain equation.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.cm",
            description="Membrane capacitance [F/m²] for the 1D Purkinje monodomain equation.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.vm1DRest",
            description=(
                "Resting transmembrane potential [V] used to initialise the 1D Purkinje field. "
                "Default: -0.084 V."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.outputVariables.export",
            description=(
                "Word list of variables written to postProcessing/purkinjeNetwork.dat. "
                "Valid tokens: Vm, Iion, Icoupling, IcouplingSource, IcouplingCurrent. "
                "Default: (Vm Icoupling)."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word_list",
            ui_control="token_list",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeNetworkModelCoeffs.outputVariables.debug",
            description=(
                "Word list of variables printed to terminal every 10 time steps. "
                "Valid tokens: Vm, Iion, Icoupling, IcouplingSource, IcouplingCurrent. "
                "Default: empty."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word_list",
            ui_control="token_list",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
            description=(
                "PVJ coupling model selector within the domainCouplings sub-dictionary. "
                "Used to couple a Purkinje network domain to the 3D myocardium."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/electroDomainCoupler.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("reactionDiffusionPvjCoupler valid only with monodomainSolver+monodomain1DSolver; eikonalPvjCoupler valid only with eikonalSolver+eikonalSolver (1D).",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj",
            description=(
                "PVJ coupling resistance [Ω·m²]. Required for "
                "reactionDiffusionPvjCoupler; unused by eikonalPvjCoupler."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjRadius",
            description=(
                "Sphere radius [m] for current distribution at PVJ junctions. "
                "Default: 0.5e-3."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
            description=(
                "PVJ coupling direction. "
                "unidirectional: 1-way, Purkinje drives tissue only. "
                "bidirectional: physiological retrograde conduction allowed for "
                "reactionDiffusionPvjCoupler. eikonalPvjCoupler currently accepts only "
                "unidirectional; bidirectional reports 'in development'."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
                "src/electroModels/electroCouplers/pvjCoupler/eikonal/eikonalPvjCoupler.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("unidirectional", "bidirectional"),
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.conductionNetworkDomain",
            description=(
                "Name of the conductionNetworkDomains entry this coupling is linked to. "
                "Required for every conduction-system coupling entry."
            ),
            source_refs=(
                "src/electroModels/core/system/electrophysicsSystemBuilder.C",
            ),
            value_kind="word",
            ui_control="text",
            dynamic_path=True,
            required=True,
            constraints=("Must match a key in conductionNetworkDomains.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme",
            description=(
                "Electrophysics time-advance scheme selector. "
                "staggeredElectrophysicsAdvanceScheme: explicit staggered, safe for unidirectional coupling. "
                "pimpleStaggeredElectrophysicsAdvanceScheme: PIMPLE-iterated, required for bidirectional coupling."
            ),
            source_refs=(
                "src/electroModels/core/advanceSchemes/electrophysicsAdvanceScheme.C",
                "src/electroModels/core/advanceSchemes/staggered/staggeredElectrophysicsAdvanceScheme.H",
                "src/electroModels/core/advanceSchemes/pimpleStaggered/pimpleStaggeredElectrophysicsAdvanceScheme.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=(
                "staggeredElectrophysicsAdvanceScheme",
                "pimpleStaggeredElectrophysicsAdvanceScheme",
            ),
            required=False,
        ),
    ),
    "active_tension": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.activeTensionModel",
            description="Active-tension model selector.",
            source_refs=(
                "src/activeTensionModels/activeTensionModel/activeTensionModel.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("GoktepeKuhl", "NashPanfilov"),
            required=False,
            constraints=("Only applicable when electro-mechanical coupling is configured.",),
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
            required=False,
            constraints=("Only applicable when activeTensionModel is configured.",),
        ),
    ),
    "domain_couplings": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
            description=(
                "Coupling model selector for domain-to-domain interactions "
                "(e.g., Purkinje-to-myocardium). Used only when multiple domains are present."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/electroDomainCoupler.C",
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("reactionDiffusionPvjCoupler valid only with monodomainSolver+monodomain1DSolver; eikonalPvjCoupler valid only with eikonalSolver+eikonalSolver (1D).",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.conductionNetworkDomain",
            description=(
                "REQUIRED: Name of the conductionNetworkDomains entry this coupling targets. "
                "Must match a key in conductionNetworkDomains."
            ),
            source_refs=(
                "src/electroModels/core/system/electrophysicsSystemBuilder.C",
            ),
            value_kind="word",
            ui_control="text",
            dynamic_path=True,
            required=True,
            constraints=("Must match a key in conductionNetworkDomains.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj",
            description=(
                "Junction resistance [Ω·m²] per unit surface area. "
                "Used by reactionDiffusionPvjCoupler for monodomain coupling. "
                "Example: 500.0 Ω·m²."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjRadius",
            description=(
                "Sphere radius [m] around each PVJ used to identify 3D cells for coupling. "
                "Default: 0.5e-3 m."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.C",
            ),
            value_kind="scalar",
            ui_control="number",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
            description=(
                "Coupling direction: 'unidirectional' (Purkinje→myocardium only) "
                "or 'bidirectional' (both ways). Default: unidirectional."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.H",
            ),
            value_kind="enum",
            ui_control="select",
            enum_values=("unidirectional", "bidirectional"),
            dynamic_path=True,
            required=False,
        ),
    ),
}


def all_documented_driver_paths() -> tuple[str, ...]:
    paths = [entry.driver_path for entry in PHYSICS_PROPERTY_ENTRIES]
    for entries in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        paths.extend(entry.driver_path for entry in entries)
    return tuple(dict.fromkeys(paths))
