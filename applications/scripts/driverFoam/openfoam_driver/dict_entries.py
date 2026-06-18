#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     dict_entries
#
# Description
#     Defines schema contracts for OpenFOAM dictionary parameters.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final, Literal
from .ionic_model_catalog import BATCHED_MODELS

# Ionic models that implement transmural tissue heterogeneity
# (configureIonicHeterogeneity, endo/M/epi blend) on CPU and/or GPU.
HETEROGENEITY_MODELS: tuple[str, ...] = (
    "BuenoOrovio", "TNNP", "TWorld", "ToRORd_dynCl",
    "BuenoOroviocompactBatched", "TNNPcompactBatched",
    "TWorldcompactBatched", "ToRORd_dynClcompactBatched",
)

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
    # Plan §5 — Structured constraints (P5a, foundation; migration in P5b).
    # Until each entry's prose ``constraints`` is migrated to one or more
    # of the structured fields below, the validator falls back to the
    # English form. All four fields default to empty, so adding an entry
    # without filling them is the additive backward-compatible case.
    #
    # Each ``{key: value}`` pair encodes a value predicate: the entry's
    # applicability/forbiddenness/requiredness is gated on ``context[key]``
    # equalling ``value`` (or appearing in the tuple when ``value`` is a
    # tuple). Block-presence predicates use virtual keys starting with
    # ``"$"`` (e.g. ``"$ecgDomains_present"``).
    applicable_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    forbidden_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    required_when: dict[str, str | tuple[str, ...]] = field(default_factory=dict)
    # Sibling-key mutual exclusion. Either side may declare the relation;
    # the validator treats it as symmetric.
    mutually_exclusive_with: tuple[str, ...] = ()


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


CONTROL_DICT_ENTRIES: Final[tuple[DictEntry, ...]] = (
    DictEntry(
        driver_path="deltaT",
        phases=frozenset({"solver"}),
        description=(
            "Simulation time step. Critical for ODE solver stability and "
            "manufactured-solution convergence tests — sweep alongside mesh "
            "refinement (number_cells) to measure temporal order. Use a large "
            "value for smoke-test runs that verify setup before committing to a "
            "fine-resolution sweep."
        ),
        source_refs=(
            "applications/solvers/cardiacFoam/cardiacFoam.C",
        ),
        value_kind="openfoam_literal",
        unit="s",
        required=True,
        typical_value="",
    ),
    DictEntry(
        driver_path="endTime",
        phases=frozenset({"solver"}),
        description=(
            "Simulation end time. Set to a small value (e.g. 1e-3) for a "
            "smoke-test run that verifies the case launches and runs at least "
            "one step without crashing, before committing to a full-length "
            "production run."
        ),
        source_refs=(
            "applications/solvers/cardiacFoam/cardiacFoam.C",
        ),
        value_kind="openfoam_literal",
        unit="s",
        required=True,
        typical_value="",
    ),
    DictEntry(
        driver_path="startTime",
        phases=frozenset({"solver"}),
        description="Simulation start time. Almost always 0 for new runs.",
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="openfoam_literal",
        unit="s",
        required=True,
        typical_value="0",
    ),
    DictEntry(
        driver_path="startFrom",
        phases=frozenset({"solver"}),
        description=(
            "Which time directory to start from. "
            "startTime uses the value of startTime; "
            "latestTime restarts from the last written time directory."
        ),
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="enum",
        enum_values=("startTime", "firstTime", "latestTime"),
        required=True,
        typical_value="startTime",
    ),
    DictEntry(
        driver_path="stopAt",
        phases=frozenset({"solver"}),
        description="Condition that halts the run.",
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="enum",
        enum_values=("endTime", "writeNow", "noWriteNow", "nextWrite"),
        required=True,
        typical_value="endTime",
    ),
    DictEntry(
        driver_path="writeControl",
        phases=frozenset({"solver"}),
        description=(
            "Trigger for writing output to disk. "
            "runTime writes every writeInterval seconds of simulation time; "
            "timeStep writes every writeInterval time steps."
        ),
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="enum",
        enum_values=("runTime", "timeStep", "clockTime", "cpuTime"),
        required=True,
        typical_value="runTime",
    ),
    DictEntry(
        driver_path="writeInterval",
        phases=frozenset({"solver"}),
        description=(
            "Output writing frequency in units of writeControl. "
            "When writeControl=runTime this is seconds of simulation time. "
            "Typical cardiac simulations write every 5 ms."
        ),
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="openfoam_literal",
        unit="s (when writeControl=runTime)",
        required=True,
        typical_value="5e-3",
    ),
    DictEntry(
        driver_path="writeFormat",
        phases=frozenset({"solver"}),
        description="Binary or ASCII output format. ASCII is human-readable; binary is faster and smaller.",
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="enum",
        enum_values=("ascii", "binary"),
        required=True,
        typical_value="ascii",
    ),
    DictEntry(
        driver_path="purgeWrite",
        phases=frozenset({"solver"}),
        description=(
            "Number of output time directories to keep on disk (0 = keep all). "
            "Use 2-3 when disk space is limited on long convergence sweeps."
        ),
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
        value_kind="openfoam_literal",
        required=True,
        typical_value="0",
    ),
)


ELECTRO_PROPERTY_ENTRY_GROUPS: Final[dict[str, tuple[DictEntry, ...]]] = {
    "top_level": (
        DictEntry(
            driver_path="cellZone",
            phases=frozenset({"anatomy"}),
            description="Optional cell zone name. Restricts the myocardium domain and equations to a specific subset of the mesh cells.",
            source_refs=("src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C",),
            value_kind="word",
            required=False,
            typical_value="heart",
        ),
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
            description=(
                "Diffusion step time-discretisation for monodomain and bidomain solvers. "
                "'explicit' uses an operator-split forward-Euler diffusion step; "
                "'implicit' uses a pimple-controlled implicit solve. "
                "Not applicable to eikonalSolver (steady BVP) or singleCellSolver (no diffusion)."
            ),
            source_refs=(
                "src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C",
            ),
            value_kind="enum",
            enum_values=("implicit", "explicit"),
            typical_value="explicit",
            required=False,
            constraints=(
                "Only applicable when myocardiumSolver is monodomainSolver or bidomainSolver.",
            ),
            required_when={"myocardiumSolver": ("monodomainSolver", "bidomainSolver")},
            forbidden_when={"myocardiumSolver": ("eikonalSolver", "singleCellSolver")},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activationThreshold",
            phases=frozenset({"solver"}),
            description=(
                "Vm threshold (volts) used to detect cell activation onset. "
                "A cell is marked as activated when Vm crosses this value upward. "
                "Default 0.0 V (rest potential for most models)."
            ),
            source_refs=(
                "src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.C",
            ),
            value_kind="scalar",
            typical_value="0.0",
            required=False,
            applicable_when={"myocardiumSolver": ("monodomainSolver", "bidomainSolver")},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicModel",
            phases=frozenset({"physics"}),
            description="Ionic cell model selector.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
            ),
            value_kind="enum",
            enum_values=tuple(
                [
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
                    "bathBidomainFDAManufactured",
                    "bidomainFDAManufactured",
                    "monodomainFDAManufactured",
                ]
                + BATCHED_MODELS
            ),
            required=True,
            constraints=("Not applicable when myocardiumSolver=eikonalSolver.",),
            forbidden_when={"myocardiumSolver": "eikonalSolver"},
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
            # Two independent forbidden predicates; each fires independently when
            # its own condition matches. The manufactured-model set covers all
            # three FDA models (mono, bi, bath-bi).
            applicable_when={
                "myocardiumSolver": ("monodomainSolver", "bidomainSolver", "singleCellSolver"),
                "ionicModel": tuple(
                    [
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
                    ]
                    + BATCHED_MODELS
                ),
            },
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.sex",
            phases=frozenset({"physics"}),
            description="Optional biological sex selector for ionic models that expose sex-specific variants.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/ionicModels/ionicModel/ionicSelector.C",
                "src/ionicModels/TWorld/TWorld.C",
            ),
            value_kind="enum",
            enum_values=("neutral", "male", "female"),
            required=False,
            applicable_when={"ionicModel": ("TWorld", "TWorldBatched")},
            constraints=("Only applicable for ionic models whose supportedSexTypes() includes the selected value.",),
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
            notes="staggeredElectrophysicsAdvanceScheme: weakly coupled, fast, stable for unidirectional. pimpleStaggeredElectrophysicsAdvanceScheme: strongly coupled with PIMPLE iteration, stable for bidirectional coupling (requires solutionAlgorithm=implicit in monodomainSolverCoeffs/bidomainSolverCoeffs).",
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
            applicable_when={"ionicModel": (
                "monodomainFDAManufactured",
                "bidomainFDAManufactured",
                "bathBidomainFDAManufactured",
            )},
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
            driver_path="$ELECTRO_MODEL_COEFFS.writeFrequency",
            phases=frozenset({"solver"}),
            description=(
                "Output period (seconds) for single-cell trace writing. "
                "When set, a row is emitted only when a writeFrequency boundary "
                "is crossed between the previous and current time. "
                "Omit (or set to 0) to write every step."
            ),
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
    "ionic_heterogeneity": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.field",
            phases=frozenset({"physics"}),
            description=(
                "Name of the transmural-distance field (0=endo, 1=epi) that "
                "drives the heterogeneity blend. Read by the myocardium domain "
                "when an ionicHeterogeneity block is present."
            ),
            source_refs=(
                "src/electroModels/electroDomains/myocardiumDomain/myocardiumDomainInterface.C",
            ),
            value_kind="word",
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
            examples=("t",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mode",
            phases=frozenset({"physics"}),
            description="Heterogeneity application mode.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="enum",
            enum_values=("transmuralBands",),
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.endoMInterface",
            phases=frozenset({"physics"}),
            description=(
                "Transmural position of the endocardium/M-cell boundary, "
                "normalized in (0, mEpiInterface). Defaults to 0.3 in the solver."
            ),
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="scalar",
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
            constraints=("Must be > 0 and < mEpiInterface.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mEpiInterface",
            phases=frozenset({"physics"}),
            description=(
                "Transmural position of the M-cell/epicardium boundary, "
                "normalized in (endoMInterface, 1). Defaults to 0.7 in the solver."
            ),
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="scalar",
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
            constraints=("Must be > endoMInterface and < 1.",),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionWidth",
            phases=frozenset({"physics"}),
            description=(
                "Width of the smooth transition band between tissue regions. "
                "0 selects hard (sharp) transitions. Defaults to 0.1 in the solver."
            ),
            source_refs=(
                "src/ionicModels/ionicModel/ionicModel.C",
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="scalar",
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode",
            phases=frozenset({"physics"}),
            description="Transition style between tissue bands.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="enum",
            enum_values=("blend", "hard"),
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing",
            phases=frozenset({"physics"}),
            description="Smoothing curve applied within transition bands.",
            source_refs=(
                "src/ionicModels/ionicModel/ionicHeterogeneity.C",
            ),
            value_kind="enum",
            enum_values=("smoothstep",),
            required=False,
            applicable_when={"ionicModel": HETEROGENEITY_MODELS},
        ),
    ),
    "ionic_constant_overrides": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.scale.<AC_name>",
            phases=frozenset({"physics"}),
            description=(
                "Scale an ionic model constant by a multiplicative factor. "
                "Applied globally (every cell). "
                "Use for drug effects, channelopathies, or ischaemia: e.g. "
                "AC_g_Kr 0.5 halves IKr conductance (LQT2/hERG block). "
                "Constant names have the 'AC_' prefix and are listed in "
                "ionic_model_catalog.py under the model's 'constants' field. "
                "scale is applied before set; using both on the same constant is an error."
            ),
            source_refs=(
                "src/genericWriter/ionicModelIO.C",
                "src/ionicModels/ionicModel/ionicModel.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            notes=(
                "Common TNNP examples: AC_g_Kr (IKr/hERG), AC_g_Ks (IKs/KCNQ1), "
                "AC_g_Na (INa/SCN5A), AC_g_CaL (ICaL/CACNA1C), AC_g_K1 (IK1). "
                "ToRORd_dynCl uses the same AC_ prefix; check the catalog for exact names."
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ionicConstantOverrides.global.set.<AC_name>",
            phases=frozenset({"physics"}),
            description=(
                "Set an ionic model constant to an absolute value, overriding the "
                "hardcoded default entirely. "
                "Applied globally (every cell). "
                "Use when a precise literature value is known and scaling from the "
                "default would be unreliable. "
                "set is applied after scale; using both on the same constant is an error."
            ),
            source_refs=(
                "src/genericWriter/ionicModelIO.C",
                "src/ionicModels/ionicModel/ionicModel.C",
            ),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
        ),
    ),
    "batched_integrator": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.batchedIntegrator",
            phases=frozenset({"solver"}),
            description="Integration scheme for GPU/SOA batched models.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="enum",
            enum_values=("euler", "rushLarsen"),
            required=True,
            applicable_when={"ionicModel": tuple(BATCHED_MODELS)},
            required_when={"ionicModel": tuple(BATCHED_MODELS)},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.batchedSubsteps",
            phases=frozenset({"solver"}),
            description="Number of substeps for batched ODE integration.",
            source_refs=("src/ionicModels/ionicModel/ionicModel.C",),
            value_kind="integer",
            required=True,
            typical_value="50",
            applicable_when={"ionicModel": tuple(BATCHED_MODELS)},
            required_when={"ionicModel": tuple(BATCHED_MODELS)},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.batchedParallelCells",
            phases=frozenset({"solver"}),
            description="Enable multithreaded cell loop execution for batched models.",
            source_refs=("src/ionicModels/ionicModel/batchedIonicModel.H",),
            value_kind="boolean",
            required=False,
            applicable_when={"ionicModel": tuple(BATCHED_MODELS)},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.batchedParallelMinCells",
            phases=frozenset({"solver"}),
            description="Minimum cell threshold to trigger parallel execution if batchedParallelCells is true.",
            source_refs=("src/ionicModels/ionicModel/batchedIonicModel.H",),
            value_kind="integer",
            required=False,
            applicable_when={"ionicModel": tuple(BATCHED_MODELS)},
            typical_value="256",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.storeBatchedAlgebraics",
            phases=frozenset({"solver"}),
            description="Force storing volatile algebraic states in memory for export/debug.",
            source_refs=("src/ionicModels/ionicModel/batchedIonicModel.H",),
            value_kind="boolean",
            required=False,
            applicable_when={"ionicModel": tuple(BATCHED_MODELS)},
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
            typical_value="RKF45",
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
            required_when={"myocardiumSolver": "singleCellSolver"},
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
            required_when={"myocardiumSolver": "singleCellSolver"},
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
            required_when={"myocardiumSolver": "singleCellSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_amplitude",
            phases=frozenset({"stimulus"}),
            description="Stimulus amplitude for the S1/S2 protocol.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            notes=(
                "Magnitude is model-class dependent. Full ionic models "
                "(TNNP, ORd, Grandi, Courtemanche, Fabbri, ...) use values "
                "around the typical_value below in pA. Phenomenological "
                "models (AlievPanfilov, BuenoOrovio) use dimensionless "
                "scaled values around 0.4-1.0; see "
                "ionic_model_catalog.IonicModelEntry.model_type to detect "
                "this case before using typical_value verbatim."
            ),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            typical_value="60",
            required_when={"myocardiumSolver": "singleCellSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim1",
            phases=frozenset({"stimulus"}),
            description="Number of S1 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            required_when={"myocardiumSolver": "singleCellSolver"},
            typical_value="3",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.stim_period_S2",
            phases=frozenset({"stimulus"}),
            description="S2 pacing cycle length.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            required_when={"myocardiumSolver": "singleCellSolver"},
            typical_value="0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.singleCellStimulus.nstim2",
            phases=frozenset({"stimulus"}),
            description="Number of S2 pulses.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="integer",
            required=True,
            constraints=("Required when myocardiumSolver=singleCellSolver.",),
            required_when={"myocardiumSolver": "singleCellSolver"},
            typical_value="0",
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
            required_when={"myocardiumSolver": ("monodomainSolver", "eikonalSolver")},
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
            required_when={"myocardiumSolver": ("monodomainSolver", "bidomainSolver")},
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
            required_when={"myocardiumSolver": ("monodomainSolver", "bidomainSolver")},
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
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMinList",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMax",
            phases=frozenset({"stimulus"}),
            description="Single-box maximum corner for PDE stimulus application.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMaxList.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMaxList",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMinList",
            phases=frozenset({"stimulus"}),
            description="Per-box minimum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMin.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMin",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMaxList",
            phases=frozenset({"stimulus"}),
            description="Per-box maximum corners for multi-region PDE stimulation.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="vector3_list",
            required=False,
            constraints=("Mutually exclusive with stimulusLocationMax.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusLocationMax",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTime",
            phases=frozenset({"stimulus"}),
            description="Shared start time for all PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTimeList.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTimeList",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTimeList",
            phases=frozenset({"stimulus"}),
            description="Per-box start times for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusStartTime.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusStartTime",
            ),
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
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDurationList",
            phases=frozenset({"stimulus"}),
            description="Per-box durations for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusDuration.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusDuration",
            ),
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
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensityList",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensityList",
            phases=frozenset({"stimulus"}),
            description="Per-box current densities for PDE stimulus boxes.",
            source_refs=("src/genericWriter/stimulusIO.C",),
            value_kind="scalar_list",
            required=False,
            constraints=("Mutually exclusive with stimulusIntensity.",),
            mutually_exclusive_with=(
                "$ELECTRO_MODEL_COEFFS.externalStimulus.stimulusIntensity",
            ),
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.type",
            phases=frozenset({"solver"}),  # TODO: revisit phases
            description="Optional myocardium-side verification hook selector.",
            source_refs=(
                "src/electroModels/core/verificationModels/electroVerificationModel.C",
                "src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.H",
                "src/verificationModels/bidomainVerification/manufacturedFDABidomainVerifier.H",
                "src/verificationModels/bathBidomainVerification/manufacturedFDABathBidomainVerifier.H",
            ),
            value_kind="enum",
            enum_values=(
                "manufacturedFDAMonodomainVerifier",
                "manufacturedFDABidomainVerifier",
                "manufacturedFDABathBidomainVerifier",
            ),
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.enforceExactFields",
            phases=frozenset({"solver"}),
            description="Forces the solver fields to strictly match the exact analytical solution every step.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.H",),
            value_kind="boolean",
            required=False,
            typical_value="false",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.initializeFields",
            phases=frozenset({"solver"}),
            description="Forces the solver fields to match the exact analytical solution at the initial time.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.H",),
            value_kind="boolean",
            required=False,
            typical_value="true",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.gamma",
            phases=frozenset({"physics"}),
            description="Spatial scaling parameter for the manufactured analytical solution.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",),
            value_kind="scalar",
            required=False,
            typical_value="1.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.alpha",
            phases=frozenset({"physics"}),
            description="Temporal scaling parameter for the manufactured analytical solution.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",),
            value_kind="scalar",
            required=False,
            typical_value="0.01",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.k",
            phases=frozenset({"physics"}),
            description="Wave vector parameter for the manufactured analytical solution.",
            source_refs=("src/verificationModels/monodomainVerification/manufacturedFDAMonodomainVerifier.C",),
            value_kind="scalar",
            required=False,
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.verificationModel.groundElectrode",
            phases=frozenset({"physics"}),
            description="Analytically pins a ground Dirichlet reference for bath-bidomain manufactured models.",
            source_refs=("src/verificationModels/bathBidomainVerification/manufacturedFDABathBidomainVerifier.C",),
            value_kind="boolean",
            required=False,
            typical_value="true",
        ),
    ),
    "bath_potential_domain": (
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathCellZones",
            phases=frozenset({"anatomy", "physics"}),
            description="Cell zones treated as conductive bath, torso, or organ tissue in the unified extracellular-potential domain.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word_list",
            required=True,
            constraints=("Required when bathPotentialDomain is configured.",),
            applicable_when={"$bathPotentialDomain_configured": True},
            required_when={"$bathPotentialDomain_configured": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.heartCellZone",
            phases=frozenset({"anatomy", "physics"}),
            description="Cell zone treated as myocardium inside the unified extracellular-potential domain.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word",
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
            typical_value="myocardium",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.bathConductivityField",
            phases=frozenset({"physics"}),
            description="Volume field name containing bath/organ conductivity values.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="word",
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
            typical_value="bodyAndOrgansConductivity",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.phiEReferenceValue",
            phases=frozenset({"physics"}),
            description="Extracellular-potential reference value used when a reference cell is applied.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
            unit="V",
            typical_value="0.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.reportSetup",
            phases=frozenset({"solver"}),
            description="Switch that logs heart/bath mesh wiring and boundary setup.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="boolean",
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
            typical_value="false",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.phiERefPoint",
            phases=frozenset({"physics"}),
            description="Point used to locate the reference cell for pure-Neumann extracellular-potential boundaries.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="vector3",
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.groundPatches.<patch>",
            phases=frozenset({"physics"}),
            description="Dirichlet ground-patch value map for extracellular potential.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.bathPotentialDomain.surfaceCurrentPatches.<patch>",
            phases=frozenset({"physics", "stimulus"}),
            description="Neumann surface-current patch map; values are scalar surface-current values.",
            source_refs=("src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$bathPotentialDomain_configured": True},
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
            required_when={"myocardiumSolver": "eikonalSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMin",
            phases=frozenset({"stimulus"}),
            description="Minimum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
            required_when={"myocardiumSolver": "eikonalSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.stimulusLocationMax",
            phases=frozenset({"stimulus"}),
            description="Maximum corner of the eikonal stimulus bounding box.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
            required_when={"myocardiumSolver": "eikonalSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.c0",
            phases=frozenset({"physics"}),
            description="Wave-speed parameter for the eikonal formulation.",
            source_refs=("src/electroModels/myocardiumModels/eikonalSolver/eikonalSolver.C",),
            value_kind="scalar",
            required=True,
            constraints=("Required when myocardiumSolver=eikonalSolver.",),
            required_when={"myocardiumSolver": "eikonalSolver"},
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
            # Virtual key: set by the driver when ecgDomains sub-dict is present.
            applicable_when={"$ecgDomains_present": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver",
            phases=frozenset({"physics"}),
            description=(
                "ECG solver selector within the ecgDomains sub-dictionary. "
                "pseudoECG: dipole approximation, works with any PDE solver. "
                "torsoECG: reads phiE from a bath/torso domain, requires bidomainSolver + bathPotentialDomain. "
                "eikonalECG: fast activation-time surrogate using precomputed tissue templates, requires eikonalSolver; "
                "a 'sampling' sub-dictionary (start, end, deltaT) is mandatory."
            ),
            source_refs=(
                "src/electroModels/electroDomains/ecgDomain/ecgSolver.C",
                "src/electroModels/ecgModels/eikonalECG/eikonalECG.C",
            ),
            value_kind="enum",
            enum_values=("pseudoECG", "torsoECG", "eikonalECG"),
            dynamic_path=True,
            required=False,
            constraints=(
                "Only applicable when ecgDomains block is present in electroProperties. "
                "torsoECG requires bathPotentialDomain with bidomainSolver. "
                "eikonalECG requires eikonalSolver and a mandatory 'sampling' sub-dict.",
            ),
            applicable_when={"$ecgDomains_present": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sigmaExtracellular",
            phases=frozenset({"physics"}),
            description="Extracellular conductivity fallback value used by pseudoECG solver if a 3D tensor field is not available.",
            source_refs=("src/electroModels/ecgModels/pseudoECGSolver/pseudoECGSolver.H",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": "pseudoECG"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.torsoSurface",
            phases=frozenset({"anatomy"}),
            description="Path to an STL file representing the torso boundary surface.",
            source_refs=("src/electroModels/ecgModels/torsoECG/torsoECG.H",),
            value_kind="string",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": "torsoECG"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sampling.start",
            phases=frozenset({"physics"}),
            description=(
                "Start of the Vm-response interpolation window [s] for eikonalECG. "
                "Usually 0. The tissue template lookup begins at this time."
            ),
            source_refs=("src/electroModels/ecgModels/eikonalECG/eikonalECG.C",),
            value_kind="openfoam_literal",
            unit="s",
            required=True,
            typical_value="0",
            dynamic_path=True,
            applicable_when={"$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": "eikonalECG"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sampling.end",
            phases=frozenset({"physics"}),
            description=(
                "End of the Vm-response interpolation window [s] for eikonalECG. "
                "Must be >= action potential duration. "
                "Human ventricular AP: ~0.3–0.5 s. "
                "Must not exceed tissue template duration (~1.0 s). "
                "Setting this shorter than APD truncates repolarisation in the ECG."
            ),
            source_refs=("src/electroModels/ecgModels/eikonalECG/eikonalECG.C",),
            value_kind="openfoam_literal",
            unit="s",
            required=True,
            typical_value="0.5",
            dynamic_path=True,
            applicable_when={"$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": "eikonalECG"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.sampling.deltaT",
            phases=frozenset({"physics"}),
            description=(
                "Time step of the Vm-response interpolation [s] for eikonalECG. "
                "Controls ECG output resolution; independent of the solver deltaT. "
                "Typical range: 0.001–0.005 s. Finer than 0.001 s has no benefit "
                "as the tissue templates are sampled at 0.1 ms resolution."
            ),
            source_refs=("src/electroModels/ecgModels/eikonalECG/eikonalECG.C",),
            value_kind="openfoam_literal",
            unit="s",
            required=True,
            typical_value="0.005",
            dynamic_path=True,
            applicable_when={"$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": "eikonalECG"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.reportElectrodeLookup",
            phases=frozenset({"solver"}),
            description="Switch controlling torsoECG electrode-to-cell lookup logging.",
            source_refs=("src/electroModels/ecgModels/torsoECG/torsoECG.C",),
            value_kind="boolean",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgSolver=torsoECG.",),
            typical_value="true",
            # ecgSolver lives under ecgDomains.<name>.ecgSolver — a dynamic path whose
            # slot key contains the resolved <name> segment, which is not known
            # statically. This constraint resists structured encoding; prose is the
            # authoritative form.  The ecgDomains presence guard is a reasonable
            # approximation since reportElectrodeLookup only appears inside ECG blocks.
            applicable_when={"$ecgDomains_present": True},
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
            applicable_when={"$ecgDomains_present": True},
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
            applicable_when={"$ecgDomains_present": True},
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
            applicable_when={"$ecgDomains_present": True},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.checkQuadratureOrders",
            phases=frozenset({"physics"}),
            description=(
                "List of quadrature orders used to compare manufactured pseudo-ECG "
                "reference convergence. Default when omitted: a single-element list [6]."
            ),
            source_refs=("src/verificationModels/ecgVerification/pseudoECGManufacturedVerifier.C",),
            value_kind="label_list",
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when ecgDomains block is present in electroProperties.",),
            applicable_when={"$ecgDomains_present": True},
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
            applicable_when={"$ecgDomains_present": True},
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
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler", "eikonalMonodomainPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when an ECG domain declares a coupling block.",),
            # Virtual key: set when any ecgDomains.<name>.coupling block is populated.
            # Approximated here as ecgDomains_present since the coupling block
            # only appears inside ECG configurations.
            applicable_when={"$ecgDomains_present": True},
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
            required_when={"myocardiumSolver": "bidomainSolver"},
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
            required_when={"myocardiumSolver": "bidomainSolver"},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiERefPoint",
            phases=frozenset({"physics"}),
            description=(
                "Point [m] used to locate the cell that pins the extracellular "
                "potential reference. The bidomain φE equation has pure Neumann "
                "boundary conditions, leaving φE determined only up to an "
                "additive constant — one cell must be clamped to break the "
                "indeterminacy. Monodomain solves only Vm with mixed BCs and "
                "does not need a reference point."
            ),
            notes=(
                "Applies to plain bidomain and bath-bidomain alike — both run "
                "on bidomainSolver, only the extracellular domain (heart-only "
                "vs heart+bath) differs. If the supplied point falls outside "
                "the mesh the solver falls back to a default cell."
            ),
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="vector3",
            required=True,
            constraints=("Required for bidomainSolver.",),
            required_when={"myocardiumSolver": "bidomainSolver"},
            typical_value="(0 0 0)",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.phiEReferenceValue",
            phases=frozenset({"physics"}),
            description=(
                "Value of the extracellular potential [V] at the reference "
                "cell located by phiERefPoint. Optional — defaults to 0 V "
                "when omitted, which is the standard choice for ungrounded "
                "bidomain solves."
            ),
            source_refs=("src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.C",),
            value_kind="scalar",
            required=False,
            unit="V",
            typical_value="0",
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
            enum_values=("purkinjeGraphModel",),
            dynamic_path=True,
            required=False,
            constraints=("Only applicable when conductionNetworkDomains block is present.",),
            # Virtual key: set by the driver when a conductionNetworkDomains sub-dict
            # is present in the run config. Dynamic-path entries (<name>) cannot be
            # evaluated against literal slot keys in a flat context.
            applicable_when={"$conductionNetworkDomains_present": True},
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
            # Virtual key: set when conductionNetworkDomains block is present.
            # Dynamic path prevents static slot-key evaluation.
            required_when={"$conductionNetworkDomains_present": True},
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
            enum_values=("monodomain1DSolver", "eikonalSolver1D", "restitutionEikonalSolver1D"),
            dynamic_path=True,
            required=False,
            constraints=("monodomain1DSolver valid only with monodomainSolver myocardium; eikonalSolver1D valid only with eikonalSolver myocardium.",),
            # This is a pairing constraint (own value ↔ sibling myocardiumSolver value).
            # Neither forbidden_when nor required_when can express "iff own value == X
            # then sibling == Y" without a cross-product of own-value and sibling-value.
            # Prose-only; no structured form added.
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.useEdgeConductance",
            phases=frozenset({"physics"}),
            description="If true, uses local edge conductances to scale conduction velocity in the eikonal solver.",
            source_refs=("src/electroModels/conductionSystemModels/eikonalSolver1D/restitutionEikonalSolver1D.C",),
            value_kind="boolean",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.conductionSystemSolver": ("restitutionEikonalSolver1D",)},
            typical_value="true",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.referenceConductance",
            phases=frozenset({"physics"}),
            description="Reference conductance value used to normalize the edge conductances when scaling velocity.",
            source_refs=("src/electroModels/conductionSystemModels/eikonalSolver1D/restitutionEikonalSolver1D.C",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.conductionSystemSolver": ("restitutionEikonalSolver1D",)},
            typical_value="1.0",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.apdNominal",
            phases=frozenset({"physics"}),
            description="Nominal Action Potential Duration [ms] for multi-beat restitution dynamics.",
            source_refs=("src/electroModels/conductionSystemModels/eikonalSolver1D/restitutionEikonalSolver1D.H",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.conductionSystemSolver": ("restitutionEikonalSolver1D",)},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.escapeInterval",
            phases=frozenset({"physics"}),
            description="Funny current escape interval [ms] dictating spontaneous firing in absence of stimulus.",
            source_refs=("src/electroModels/conductionSystemModels/eikonalSolver1D/restitutionEikonalSolver1D.H",),
            value_kind="scalar",
            dynamic_path=True,
            required=False,
            applicable_when={"$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>.purkinjeGraphModelCoeffs.conductionSystemSolver": ("restitutionEikonalSolver1D",)},
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
            enum_values=tuple(
                [
                    "AlievPanfilov",
                    "BuenoOrovio",
                    "Courtemanche",
                    "Fabbri",
                    "Gaur",
                    "Grandi",
                    "Stewart",
                    "TNNP",
                    "ToRORd_dynCl",
                    "Trovato",
                    "monodomainFDAManufactured",
                    "bidomainFDAManufactured",
                ]
                + BATCHED_MODELS
            ),
            notes="Stewart is the canonical human Purkinje model. monodomainFDAManufactured/bidomainFDAManufactured are for verification only.",
            dynamic_path=True,
            required=False,
            constraints=("Required when conductionSystemSolver=monodomain1DSolver.",),
            # conductionSystemSolver lives at a dynamic path (sibling key with <name>)
            # so its slot key in a flat context contains literal "<name>" — it will
            # never match a real run's resolved key. Prose-only for this predicate;
            # the virtual-key guard is a reasonable documentation placeholder.
            required_when={"$conductionNetworkDomains_present": True},
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
                "are a hardcoded set: Vm, Iion, activationTime, "
                "IcouplingSource, IcouplingCurrent. Unknown tokens are "
                "silently dropped without a warning — see the if/else "
                "ladder in conductionSystemDomain.C. Default: "
                "(Vm IcouplingSource)."
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
                "activationTime, IcouplingSource, IcouplingCurrent. "
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
            enum_values=(
                "GoktepeKuhl", "NashPanfilov", "LandNiederer",
                "GoktepeKuhlBatched", "NashPanfilovBatched", "LandNiedererBatched",
                "ManufacturedElectromechanics",
            ),
            required=False,
            constraints=("Only applicable when electro-mechanical coupling is configured.",),
            # Virtual key: set when physics.type=electroMechanicalModel. The top-level
            # "type" key is owned by PHYSICS_PROPERTY_ENTRIES; its slot key is "type".
            applicable_when={"type": "electroMechanicalModel"},
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
            # activeTensionModel slot key is "activeTensionModel.activeTensionModel"
            # (post-prefix strip). Checking for a non-empty value signals presence.
            applicable_when={"activeTensionModel.activeTensionModel": ("GoktepeKuhl", "NashPanfilov")},
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.activeTensionModel.TaScale",
            phases=frozenset({"physics"}),
            description="Scales the active tension magnitude computed by the model.",
            source_refs=("src/activeTensionModels/activeTensionModel/activeTensionModel.C",),
            value_kind="scalar",
            required=False,
            applicable_when={"activeTensionModel.activeTensionModel": ("GoktepeKuhl", "NashPanfilov")},
            typical_value="1e3",
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
                "src/electroModels/electroCouplers/pvjCoupler/eikonalMonodomain/eikonalMonodomainPvjCoupler.C",
            ),
            value_kind="enum",
            enum_values=("reactionDiffusionPvjCoupler", "eikonalPvjCoupler", "eikonalMonodomainPvjCoupler"),
            dynamic_path=True,
            required=False,
            constraints=("reactionDiffusionPvjCoupler valid only with monodomainSolver+monodomain1DSolver; eikonalPvjCoupler valid only with eikonalSolver+eikonalSolver (1D).",),
            # Pairing constraint (own value ↔ pair of sibling keys across two separate
            # domains). Not encodable as a single forbidden_when/required_when predicate
            # without own-value introspection. Prose-only; no structured form added.
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
            # Referential-integrity constraint (value must equal a key in a sibling
            # block). Not encodable in any of the four structured families.
            # Prose-only; no structured form added.
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
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjKernel",
            phases=frozenset({"physics"}),
            description=(
                "Interpolation kernel used when mapping 1D Purkinje junctions to 3D myocardium. "
                "Options: uniform, gaussian, linear. Default: uniform."
            ),
            source_refs=(
                "src/electroModels/electroCouplers/pvjCoupler/pvjCoupler.C",
            ),
            value_kind="enum",
            enum_values=("uniform", "gaussian", "linear"),
            dynamic_path=True,
            required=False,
            typical_value="uniform",
        ),
        DictEntry(
            driver_path="$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.debugCoupling",
            phases=frozenset({"physics"}),
            description="Outputs verbose logging about PVJ coupling weights and mappings.",
            source_refs=("src/electroModels/electroCouplers/electroDomainCoupler.C",),
            value_kind="boolean",
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
