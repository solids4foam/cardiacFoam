"""cardiacFoam-owned physicsProperties and controlDict catalog entries.

The legacy ``openfoam_driver.dict_entries`` module re-exports these values so
Plan 1 changes ownership without changing import paths or serialized catalogs.
"""

from __future__ import annotations

from typing import Final

from openfoam_driver.core.contracts.dictionary import DictEntry


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
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
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
        source_refs=("applications/solvers/cardiacFoam/cardiacFoam.C",),
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
