from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal

# Workflow phases used by run documents and catalog exports, in strict order.
# Every ``DictEntry`` may declare one or more of these in ``phases``; the
# catalog exporter fans it out to each phase bucket so multi-phase entries
# appear in every consumer that needs to see them. The *primary* phase is
# resolved at validation time as the first phase in this order that the
# entry claims.
#
# Post-completion analysis (formerly the "review" phase) lives in the
# Reports section of the workspace, not in the phase walk. See spec
# section 5a and the report_catalog manifest.
Phase = Literal["anatomy", "physics", "stimulus", "solver"]


@dataclass(frozen=True)
class DictEntry:
    driver_path: str
    description: str
    source_refs: tuple[str, ...]
    notes: str = ""
    value_kind: str = "openfoam_literal"
    enum_values: tuple[str, ...] = ()
    examples: tuple[str, ...] = ()
    dynamic_path: bool = False
    required: bool = False
    constraints: tuple[str, ...] = ()
    unit: str = ""
    typical_value: str = ""
    # Workflow phases this entry belongs to. Empty is allowed transiently
    # during migration (Task A2 classifies every entry); a coverage test
    # enforces non-empty once classification lands. Values are drawn from
    # the ``Phase`` literal.
    phases: frozenset[str] = frozenset()


PHYSICS_PROPERTY_ENTRIES: Final[tuple[DictEntry, ...]] = (
    DictEntry(
        driver_path="type",
        phases=frozenset({"physics"}),
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
        enum_values=("electroModel", "electroMechanicalModel"),
        required=True,
    ),
)


ELECTRO_PROPERTY_ENTRY_GROUPS: Final[dict[str, tuple[DictEntry, ...]]] = {
    "top_level": (
        DictEntry(
            driver_path="myocardiumSolver",
            phases=frozenset({"physics"}),
            description="Top-level myocardium solver selector. Determines the active '<solver>Coeffs' sub-dictionary in electroProperties.",
            source_refs=(
                "src/electroModels/core/electroModel.C",
                "src/electroModels/core/electrophysiologyModel/electrophysiologyModel.C",
                "src/electroModels/electroDomains/myocardiumDomain/myocardiumSolver.C",
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",
                "src/electroModels/myocardiumModels/singleCellSolver/singleCellSolver.C",
                "src/electroModels/electroDomains/myocardiumDomain/eikonalMyocardiumDomain.C",
            ),
            value_kind="enum",
            enum_values=("monodomainSolver", "bidomainSolver", "singleCellSolver", "eikonalSolver"),
            notes=(
                "Dictionary values are top-level myocardiumSolver choices. "
                "monodomainSolver and bidomainSolver dispatch through the myocardiumSolver RTST "
                "inside electrophysiologyModel; singleCellSolver is registered in the parent "
                "electroModel RTST; eikonalSolver selects the electrophysiologyModel wrapper and "
                "then builds an eikonalMyocardiumDomain, not the legacy 3D eikonalSolver class."
            ),
            required=True,
        ),
    ),
    "common_model_coeffs": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solutionAlgorithm",
            phases=frozenset({"solver"}),
            description="Electro solver time-discretisation mode.",
            source_refs=(
                "src/electroModels/core/electroModel.C",
                "src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",
            ),
            value_kind="enum",
            enum_values=("implicit", "explicit"),
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicModel",
            phases=frozenset({"physics"}),
            description="Ionic cell model selector.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
            ),
            value_kind="enum",
            enum_values=(
                "AlievPanfilov",
                "BuenoOrovio",
                "Courtemanche",
                "Fabbri",
                "Gaur",
                "Grandi",
                "ORd",
                "PerisYague",
                "Stewart",
                "TNNP",
                "ToRORd_dynCl",
                "Trovato",
                "TWorld",
                "bathBidomainFDAManufactured",
                "bidomainFDAManufactured",
                "monodomainFDAManufactured",
            ),
            required=True,
            constraints=("Not applicable when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.tissue",
            phases=frozenset({"physics"}),
            description="Tissue selector for ionic-model specialisation.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicSelector.C",
            ),
            value_kind="enum",
            enum_values=("epicardialCells", "mCells", "endocardialCells", "myocyte"),
            required=True,
            constraints=("Not applicable when myocardiumSolver=eikonalSolver or ionicModel is a manufactured model.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme",
            phases=frozenset({"solver"}),
            description="Time-advance scheme for multi-domain coupling (myocardium, Purkinje, ECG).",
            source_refs=(
                "src/electroModels/core/advanceSchemes/electrophysicsAdvanceScheme.H",
                "src/electroModels/core/advanceSchemes/staggered/staggeredElectrophysicsAdvanceScheme.C",
                "src/electroModels/core/advanceSchemes/pimpleStaggered/pimpleStaggeredElectrophysicsAdvanceScheme.C",
            ),
            notes="staggeredElectrophysicsAdvanceScheme: weakly coupled, fast, stable for unidirectional. pimpleStaggeredElectrophysicsAdvanceScheme: strongly coupled with PIMPLE iteration, stable for bidirectional coupling (requires solutionAlgorithm=implicit).",
            value_kind="enum",
            enum_values=("staggeredElectrophysicsAdvanceScheme", "pimpleStaggeredElectrophysicsAdvanceScheme"),
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.dimension",
            phases=frozenset({"physics"}),
            description="Dimensional selector used by manufactured/verification ionic models.",
            source_refs=("src/ionicModels/ionicModel/ionicSelector.C",),
            value_kind="enum",
            enum_values=("1D", "2D", "3D"),
            required=False,
            constraints=("Only applicable for manufactured ionic models (monodomainFDAManufactured, bidomainFDAManufactured, bathBidomainFDAManufactured).",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.writeAfterTime",
            phases=frozenset({"solver"}),
            description="Suppresses single-cell trace output before the given time.",
            source_refs=("src/genericWriter/ionicModelIO.C",),
            value_kind="scalar",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.utilities",
            phases=frozenset({"physics"}),  # TODO: revisit phases
            description="Utility mode flag for ionic-model helper applications.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="boolean",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initSampleCell",
            phases=frozenset({"physics"}),  # TODO: revisit phases
            description="Integration-point index used for sampled single-cell style output.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="integer",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.export",
            phases=frozenset({"solver"}),
            description=(
                "Ionic variables exported to volumetric fields or trace output. "
                "Names are filtered against the active ionic model's state and "
                "algebraic-variable lists (model-dependent). 'Vm' and 'Iion' are "
                "aliased universally (resolve to the model's voltage state and "
                "total ionic current via ionicVariableCompatibility.C). Unknown "
                "names are silently dropped with a runtime warning at "
                "ionicModelIO.C:417. NOT consumed by eikonalSolver (no ionic "
                "model owned). 'activationTime' is NOT a valid ionic variable; "
                "it is a Purkinje conduction-system field."
            ),
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/genericWriter/ionicModelIO.C",
                "src/genericWriter/ionicVariableCompatibility.C",
            ),
            value_kind="word_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.ionic.debug",
            phases=frozenset({"solver"}),
            description=(
                "Ionic variables printed in debug output. Same filtering and "
                "aliasing rules as outputVariables.ionic.export."
            ),
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/genericWriter/ionicModelIO.C",
            ),
            value_kind="word_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.export",
            phases=frozenset({"solver"}),
            description="Active-tension variables exported to fields or traces.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.outputVariables.activeTension.debug",
            phases=frozenset({"solver"}),
            description="Active-tension variables printed in debug output.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.H",),
            value_kind="word_list",
            required=False,
        ),
    ),
    "ode_solver_passthrough": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.solver",
            phases=frozenset({"solver"}),
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
            enum_values=("RKF45",),
            required=True,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.initialODEStep",
            phases=frozenset({"solver"}),
            description="Initial ODE step-size hint passed through to the selected ODE solver.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            required=True,
            unit="s",
            typical_value="1e-5",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.maxSteps",
            phases=frozenset({"solver"}),
            description="Maximum internal ODE steps allowed per macro time step.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="integer",
            required=True,
            typical_value="1000",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.absTol",
            phases=frozenset({"solver"}),
            description="Absolute tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.relTol",
            phases=frozenset({"solver"}),
            description="Relative tolerance passed through to adaptive ODE solvers.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.H",),
            notes="Pass-through key; commonly used in repository tutorials.",
            value_kind="scalar",
            required=False,
        ),
    ),
    "single_cell_stimulus": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_start",
            phases=frozenset({"stimulus"}),
            description="Start time for the S1/S2 stimulus protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            unit="s",
            typical_value="0.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S1",
            phases=frozenset({"stimulus"}),
            description="S1 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            unit="s",
            typical_value="1.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_duration",
            phases=frozenset({"stimulus"}),
            description="Pulse duration for both S1 and S2 stimuli.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            unit="s",
            typical_value="1.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
            phases=frozenset({"stimulus"}),
            description="Stimulus amplitude for the S1/S2 protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            typical_value="0.4",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim1",
            phases=frozenset({"stimulus"}),
            description="Number of S1 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S2",
            phases=frozenset({"stimulus"}),
            description="S2 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim2",
            phases=frozenset({"stimulus"}),
            description="Number of S2 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
        ),
    ),
    "monodomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivity",
            phases=frozenset({"physics"}),
            description="Monodomain conductivity tensor.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            required=True,
            constraints=("Required for monodomainSolver and eikonalSolver; not used by singleCellSolver.",),
            unit="S/m",
            typical_value="0.17",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.chi",
            phases=frozenset({"physics"}),
            description="Surface-to-volume ratio.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            required=True,
            constraints=("Required for monodomainSolver and bidomainSolver; not used by singleCellSolver.",),
            unit="1/m",
            typical_value="140000",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.cm",
            phases=frozenset({"physics"}),
            description="Membrane capacitance.",
            source_refs=(
                "src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",
                "src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",
            ),
            value_kind="scalar",
            required=True,
            constraints=("Required for monodomainSolver and bidomainSolver; not used by singleCellSolver.",),
            unit="F/m²",
            typical_value="0.01",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.infoFrequency",
            phases=frozenset({"solver"}),
            description="Logging cadence for monodomain explicit stepping.",
            source_refs=("src/electroModels/myocardiumModels/monodomainSolver/monodomainSolver.C",),
            value_kind="integer",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMin",
            phases=frozenset({"stimulus"}),
            description="Single-box minimum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMinList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMax",
            phases=frozenset({"stimulus"}),
            description="Single-box maximum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMaxList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMinList",
            phases=frozenset({"stimulus"}),
            description="Per-box minimum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMin.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMaxList",
            phases=frozenset({"stimulus"}),
            description="Per-box maximum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMax.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTime",
            phases=frozenset({"stimulus"}),
            description="Shared start time for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTimeList.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTimeList",
            phases=frozenset({"stimulus"}),
            description="Per-box start times for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTime.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDuration",
            phases=frozenset({"stimulus"}),
            description="Shared duration for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            required=False,
            constraints=("Mutually exclusive with stimulusDurationList.",),
            unit="s",
            typical_value="0.002",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
            phases=frozenset({"stimulus"}),
            description="Per-box durations for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusDuration.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity",
            phases=frozenset({"stimulus"}),
            description="Shared current density for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes="DimensionedScalar-style values are best passed as a full literal string.",
            value_kind="dimensioned_scalar_literal",
            required=False,
            constraints=("Mutually exclusive with stimulusIntensityList.",),
            unit="A/m³",
            typical_value="50000",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensityList",
            phases=frozenset({"stimulus"}),
            description="Per-box current densities for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusIntensity.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.type",
            phases=frozenset({"solver"}),  # TODO: revisit phases
            description="Optional myocardium-side verification hook selector.",
            source_refs=(
                "src/verificationModels/electroVerification/electroVerificationModel.C",
                "src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.H",
                "src/verificationModels/bidomainVerification/manufacturedFDABidomainVerifier.H",
                "src/verificationModels/bidomainVerification/singleCellManufacturedFDABidomainVerifier.H",
            ),
            value_kind="enum",
            enum_values=(
                "manufacturedFDAMonodomainVerifier",
                "manufacturedFDABidomainVerifier",
                "manufacturedFDABathBidomainVerifier",
                "singleCellManufacturedFDABidomainVerifier",
            ),
            required=False,
        ),
    ),
    "potential_domain": (
        DictEntry(
            driver_path="potentialDomain.type",
            phases=frozenset({"physics"}),
            description="Top-level extracellular-potential domain selector for bidomain+bath solves.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.H",),
            value_kind="enum",
            enum_values=("extracellularPotentialDomain",),
            required=False,
            constraints=("Required when bathECGProbe or bidomainSolver uses a unified heart+bath potential domain.",),
        ),
        DictEntry(
            driver_path="potentialDomain.bathCellZones",
            phases=frozenset({"anatomy", "physics"}),
            description="Cell zones treated as conductive bath tissue in the unified extracellular-potential domain.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word_list",
            required=True,
            constraints=("Required when potentialDomain is configured.",),
        ),
        DictEntry(
            driver_path="potentialDomain.heartCellZone",
            phases=frozenset({"anatomy", "physics"}),
            description="Cell zone treated as myocardium inside the unified extracellular-potential domain.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word",
            required=False,
            typical_value="myocardium",
        ),
        DictEntry(
            driver_path="potentialDomain.bathConductivityField",
            phases=frozenset({"physics"}),
            description="Volume field name containing bath/organ conductivity values.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word",
            required=False,
            typical_value="bodyAndOrgansConductivity",
        ),
        DictEntry(
            driver_path="potentialDomain.phiEReferenceValue",
            phases=frozenset({"physics"}),
            description="Extracellular-potential reference value used when a reference cell is applied.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            required=False,
            unit="V",
            typical_value="0.0",
        ),
        DictEntry(
            driver_path="potentialDomain.reportSetup",
            phases=frozenset({"solver"}),
            description="Switch that logs heart/bath mesh wiring and boundary setup.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="boolean",
            required=False,
            typical_value="false",
        ),
        DictEntry(
            driver_path="potentialDomain.phiERefPoint",
            phases=frozenset({"physics"}),
            description="Point used to locate the reference cell for pure-Neumann extracellular-potential boundaries.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="vector3",
            required=False,
        ),
        DictEntry(
            driver_path="potentialDomain.groundPatches.<patch>",
            phases=frozenset({"physics"}),
            description="Dirichlet ground-patch value map for extracellular potential.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="potentialDomain.surfaceCurrentPatches.<patch>",
            phases=frozenset({"physics", "stimulus"}),
            description="Neumann surface-current patch map; values are scalar surface-current values.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
    ),
    "eikonal_diffusion": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.eikonalAdvectionDiffusionApproach",
            phases=frozenset({"physics"}),
            description="Toggle between the advection-diffusion and deferred-correction forms.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="boolean",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMin",
            phases=frozenset({"stimulus"}),
            description="Minimum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMax",
            phases=frozenset({"stimulus"}),
            description="Maximum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.c0",
            phases=frozenset({"physics"}),
            description="Wave-speed parameter for the eikonal formulation.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
        ),
    ),
    "ecg": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.electrodePositions.<electrode>",
            phases=frozenset({"physics"}),
            description="Shared electrode position vector [m] inherited by ECG domain entries unless overridden.",
            source_refs=(
                "src/electroModels/core/system/electrophysicsSystemBuilder.C",
                "src/electroModels/electroDomains/ecgDomain/ecgDomain.C",
            ),
            value_kind="vector3",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            phases=frozenset({"physics"}),
            description="ECG solver selector within the ecgDomains sub-dictionary.",
            source_refs=(
                "src/electroModels/electroDomains/ecgDomain/ecgSolver.C",
            ),
            value_kind="enum",
            enum_values=("bathECGProbe", "pseudoECG"),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties. bathECGProbe requires potentialDomain with bidomainSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.reportElectrodeLookup",
            phases=frozenset({"solver"}),
            description="Switch controlling bathECGProbe electrode-to-cell lookup logging.",
            source_refs=("src/electroModels/ecgModels/bathECGProbe/bathECGProbe.C",),
            value_kind="boolean",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgSolver=bathECGProbe.",),
            typical_value="true",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.enabled",
            phases=frozenset({"physics"}),
            description="Enable manufactured pseudo-ECG verification for the selected ECG domain.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="boolean",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.dimension",
            phases=frozenset({"physics"}),
            description="Dimensional selector used by the manufactured pseudo-ECG reference.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="enum",
            enum_values=("1D", "2D", "3D"),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.referenceQuadratureOrder",
            phases=frozenset({"physics"}),
            description="Quadrature order used for the manufactured pseudo-ECG reference integral.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="integer",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.checkQuadratureOrders",
            phases=frozenset({"physics"}),
            description="Additional quadrature orders used to compare manufactured pseudo-ECG reference convergence.",
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="label_list",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.electrodePositions.<electrode>",
            phases=frozenset({"physics"}),
            description="Per-electrode position vector [m] in the ECG domain block.",
            source_refs=("src/electroModels/electroDomains/ecgDomain/ecgDomain.C",),
            notes="The driver can update existing electrode entries but does not insert new electrode names.",
            value_kind="vector3",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.coupling.electroDomainCoupler",
            phases=frozenset({"physics"}),
            description="Optional ECG-domain coupling selector dispatched through electroDomainCoupler.",
            source_refs=(
                "src/electroModels/core/system/electrophysicsSystemBuilder.C",
                "src/electroModels/electroCouplers/electroDomainCoupler.C",
            ),
            value_kind="enum",
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when an ECG domain declares a coupling block.",),
        ),
    ),
    "bidomain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityIntracellular",
            phases=frozenset({"physics"}),
            description="Intracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            required=True,
            constraints=("Required for bidomainSolver.",),
            unit="S/m",
            typical_value="0.17",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductivityExtracellular",
            phases=frozenset({"physics"}),
            description="Extracellular conductivity tensor for the bidomain formulation.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            notes="Tensor entries are best overridden with a full OpenFOAM literal string.",
            value_kind="dimensioned_tensor_literal",
            required=True,
            constraints=("Required for bidomainSolver.",),
            unit="S/m",
            typical_value="0.62",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiERefPoint",
            phases=frozenset({"physics"}),
            description="Point [m] used to locate the cell that pins the extracellular potential reference.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required for bidomainSolver.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceValue",
            phases=frozenset({"physics"}),
            description="Value of the extracellular potential at the reference cell.",
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="scalar",
            required=False,
        ),
    ),
    "conduction_system": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.conductionSystemDomain",
            phases=frozenset({"anatomy", "physics"}),
            description=(
                "Conduction system domain model selector (sub-dictionary within "
                "conductionNetworkDomains). Used exclusively in Purkinje/1D pipelines."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="enum",
            enum_values=("purkinjeGraphModel", "conductionSystemDomain"),
            notes="conductionSystemDomain is a legacy compatibility alias; prefer purkinjeGraphModel for new dictionaries.",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when conductionNetworkDomains block is present.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.graphFile",
            phases=frozenset({"anatomy", "physics"}),
            description=(
                "REQUIRED: Name of the graph file in constant/ directory "
                "(e.g., 'purkinjeGraph'). The graph file must contain 'edges', 'points', "
                "'pvjNodes', 'pvjLocations' dictionaries/lists."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word",
            dynamic_path=True,
            required=True,
            constraints=("Required when a Purkinje graph is configured.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.conductionSystemSolver",
            phases=frozenset({"physics"}),
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
            enum_values=("monodomain1DSolver", "eikonalSolver1D"),
            dynamic_path=True,
            required=False,
            constraints=("monodomain1DSolver valid only with monodomainSolver myocardium; eikonalSolver1D valid only with eikonalSolver myocardium.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.ionicModel",
            phases=frozenset({"physics"}),
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
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.tissue",
            phases=frozenset({"physics"}),
            description=(
                "Tissue type for 1D Purkinje ionic model. "
                "Examples: epicardialCells, endocardialCells, mCells, myocyte."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
                "src/ionicModels/ionicModel/ionicModel.H",
            ),
            value_kind="enum",
            enum_values=("epicardialCells", "endocardialCells", "mCells", "myocyte"),
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.rootStimulus.startTime",
            phases=frozenset({"stimulus"}),
            description="Start time [s] of the root-node stimulus applied to Purkinje node 0.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.rootStimulus.duration",
            phases=frozenset({"stimulus"}),
            description="Duration [s] of the root-node stimulus pulse.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.rootStimulus.intensity",
            phases=frozenset({"stimulus"}),
            description="Amplitude [A/m³] of the root-node stimulus current applied to Purkinje node 0.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.rootStimulus.node",
            phases=frozenset({"stimulus"}),
            description="Optional graph-node index that overrides the graph-file rootNode for the applied root stimulus.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="label",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.chi",
            phases=frozenset({"physics"}),
            description="Surface-to-volume ratio [1/m] for the 1D Purkinje monodomain equation.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            unit="1/m",
            typical_value="140000",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.cm",
            phases=frozenset({"physics"}),
            description="Membrane capacitance [F/m²] for the 1D Purkinje monodomain equation.",
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            unit="F/m²",
            typical_value="0.01",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.vm1DRest",
            phases=frozenset({"physics"}),
            description=(
                "Resting transmembrane potential [V] used to initialise the 1D Purkinje field. "
                "Default: -0.084 V."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.outputVariables.export",
            phases=frozenset({"solver"}),
            description=(
                "Word list of variables written to "
                "postProcessing/purkinjeNetwork.dat. Note the FLAT layout — "
                "this dict has 'export' and 'debug' directly inside "
                "outputVariables, with NO 'ionic' sub-block (unlike the "
                "myocardium-side outputVariables.ionic.export). Valid tokens "
                "are a hardcoded set: Vm, Iion, activationTime, Icoupling, "
                "IcouplingSource, IcouplingCurrent. (Icoupling and "
                "IcouplingSource alias to the same column.) Unknown tokens "
                "are silently dropped without a warning — see the if/else "
                "ladder at conductionSystemDomain.C:253. Default: "
                "(Vm Icoupling)."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word_list",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.outputVariables.debug",
            phases=frozenset({"solver"}),
            description=(
                "Word list of variables printed to terminal every 10 time "
                "steps. Same flat layout and same valid-token set as the "
                "Purkinje outputVariables.export entry: Vm, Iion, "
                "activationTime, Icoupling, IcouplingSource, IcouplingCurrent. "
                "Default: empty."
            ),
            source_refs=(
                "src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C",
            ),
            value_kind="word_list",
            dynamic_path=True,
            required=False,
        ),
    ),
    "active_tension": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.activeTensionModel",
            phases=frozenset({"physics"}),
            description="Active-tension model selector.",
            source_refs=(
                "src/activeTensionModels/activeTensionModel/activeTensionModel.C",
            ),
            value_kind="enum",
            enum_values=("GoktepeKuhl", "NashPanfilov"),
            required=False,
            constraints=("Only applicable when electro-mechanical coupling is configured.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.couplingSignal",
            phases=frozenset({"physics"}),
            description="Coupling signal requested by the active-tension model.",
            source_refs=(
                "src/activeTensionModels/GoktepeKuhl/GoktepeKuhl.C",
                "src/activeTensionModels/NashPanfilov/NashPanfilov.C",
            ),
            value_kind="enum",
            enum_values=("Vm",),
            required=False,
            constraints=("Only applicable when activeTensionModel is configured.",),
        ),
    ),
    "domain_couplings": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler",
            phases=frozenset({"physics"}),
            description=(
                "Coupling model selector for domain-to-domain interactions "
                "(e.g., Purkinje-to-myocardium). Used only when multiple domains are present."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/electroDomainCoupler.C",
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
            ),
            value_kind="enum",
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("reactionDiffusionPvjCoupler valid only with monodomainSolver+monodomain1DSolver; eikonalPvjCoupler valid only with eikonalSolver+eikonalSolver (1D).",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.conductionNetworkDomain",
            phases=frozenset({"physics"}),
            description=(
                "REQUIRED: Name of the conductionNetworkDomains entry this coupling targets. "
                "Must match a key in conductionNetworkDomains."
            ),
            source_refs=(
                "src/electroModels/core/system/electrophysicsSystemBuilder.C",
            ),
            value_kind="word",
            dynamic_path=True,
            required=True,
            constraints=("Must match a key in conductionNetworkDomains.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj",
            phases=frozenset({"physics"}),
            description=(
                "Junction resistance [Ω·m²] per unit surface area. "
                "Used by reactionDiffusionPvjCoupler for monodomain coupling. "
                "Example: 500.0 Ω·m²."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjRadius",
            phases=frozenset({"physics"}),
            description=(
                "Sphere radius [m] around each PVJ used to identify 3D cells for coupling. "
                "Default: 0.5e-3 m."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
            phases=frozenset({"physics"}),
            description=(
                "Coupling direction: 'unidirectional' (Purkinje→myocardium only) "
                "or 'bidirectional' (both ways). Default: unidirectional."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.H",
            ),
            value_kind="enum",
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
