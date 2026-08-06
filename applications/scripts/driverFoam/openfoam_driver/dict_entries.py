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
from typing import Any, Final, Literal
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import BATCHED_MODELS

# Ionic models that implement transmural tissue heterogeneity
# (configureIonicHeterogeneity, endo/M/epi blend and/or namedRegions) on
# CPU and/or GPU.

def get_heterogeneity_models() -> tuple[str, ...]:
    from openfoam_driver.core.plugin_interface import get_active_plugin
    return get_active_plugin().get_capabilities().get("heterogeneity_models", ())

def get_electro_property_entry_groups() -> dict[str, tuple[DictEntry, ...]]:
    from openfoam_driver.core.plugin_interface import get_active_plugin
    return get_active_plugin().get_dict_groups()

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
    source_refs: tuple[str, ...] = ()
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


from dataclasses import replace

def build_group(defaults: dict[str, Any], entries: tuple[DictEntry, ...]) -> tuple[DictEntry, ...]:
    out = []
    for e in entries:
        changes = {}
        for k, v in defaults.items():
            current = getattr(e, k)
            if not current:
                changes[k] = v
            elif isinstance(current, dict) and isinstance(v, dict):
                changes[k] = {**v, **current}
        if changes:
            out.append(replace(e, **changes))
        else:
            out.append(e)
    return tuple(out)

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




def all_documented_driver_paths() -> tuple[str, ...]:
    paths = [entry.driver_path for entry in PHYSICS_PROPERTY_ENTRIES]
    for entries in get_electro_property_entry_groups().values():
        paths.extend(entry.driver_path for entry in entries)
    return tuple(dict.fromkeys(paths))
