"""Curated registry of the canonical regression-equivalence cases.

This registry is owned by driverFOAM, not auto-derived from every upstream
cardiacFoam tutorial that happens to ship a ``regression/*.reference`` file.
Upstream can add new tutorials without breaking this suite; we intentionally
expand the registry only when we want driverFOAM to take ownership of a new
equivalence case.
"""
from __future__ import annotations

from dataclasses import dataclass
@dataclass(frozen=True)
class RegressionCase:
    # Path under tutorials/, matching an Alltest-regression REGRESSION_TESTS entry.
    case_dir: str
    # Registered agent entry name, or None for cases with no agent spec.
    entry_name: str | None
    # Case-authored dicts the agent regenerates (relpaths under the case dir).
    # Used by the round-trip stability check for mapped cases.
    dicts: tuple[str, ...]
    # Representative reference-comparison file, relative to the case dir.
    reference_file: str
    regression_script: str = "regression/regressionTest.sh"
    # Whether the agent's case discovery can address the case by folder path.
    # False for layouts the agent does not recognize (e.g. electromechanical
    # cases keep electroProperties at constant/electro/, but discovery requires
    # constant/electroProperties).
    generic_addressable: bool = True

    @property
    def mapped(self) -> bool:
        return self.entry_name is not None

    @property
    def drivers(self) -> tuple[str, ...]:
        # Mapped cases are driven both ways to confirm the agent can reason
        # about them through the strict entry and the generic case-folder path.
        return ("strict", "generic") if self.mapped else ("generic",)


_ELECTRO = "constant/electroProperties"
_PHYSICS = "constant/physicsProperties"


_KNOWN_CASES: tuple[RegressionCase, ...] = (
    RegressionCase(
        "electrophysiologyProtocols/singleCell", "singleCell",
        (_ELECTRO, _PHYSICS), "regression/singleCell.reference",
    ),
    RegressionCase(
        "NiedererEtAl2011/NiedererEtAl2011verification", "niederer2012",
        (_ELECTRO, _PHYSICS), "regression/NiedererEtAl2012.reference",
    ),
    RegressionCase(
        "manufacturedSolutions/bidomain", "manufacturedBidomain",
        (_ELECTRO, _PHYSICS), "regression/bidomainManufactured.reference",
    ),
    RegressionCase(
        "manufacturedSolutions/monodomainPseudoECG", "manufacturedMonodomainPseudoECG",
        (_ELECTRO, _PHYSICS), "regression/monodomainPseudoECG.reference",
    ),
    RegressionCase(
        "manufacturedSolutions/eikonalECG", "manufacturedEikonalECG",
        (_ELECTRO, _PHYSICS), "regression/eikonalECG.reference",
    ),
    RegressionCase(
        "manufacturedSolutions/bathBidomain", "manufacturedBathBidomain",
        (_ELECTRO, _PHYSICS), "regression/bathBidomainManufactured.reference",
    ),
    RegressionCase(
        "NiedererEtAl2011/electroMechanicalNiedererEtAl2011", None,
        (), "regression/electroMechHeterogeneity.reference",
    ),
    # purkinje ships two reference files (eikonalSlab.reference,
    # purkinjeSlab.reference); reference_file records the representative one.
    RegressionCase(
        "NiedererEtAl2011/purkinjeNiedererEtAl2011", None,
        (), "regression/purkinjeSlab.reference",
    ),
    RegressionCase(
        "electrophysiologyProtocols/rotorInstability", None,
        (), "regression/rotorInstability.reference",
    ),
)
REGRESSION_CASES: tuple[RegressionCase, ...] = _KNOWN_CASES
