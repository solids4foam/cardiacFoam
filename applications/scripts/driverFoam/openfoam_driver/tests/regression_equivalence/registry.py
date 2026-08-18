"""Declarative registry of the 9 canonical cardiacFoam regressions.

Keeps the harness in lockstep with ``tutorials/Alltest-regression``. Each row
records how the driverFOAM agent should reproduce that case.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from openfoam_driver.specs.common import tutorials_root_default


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


def _discover_cases() -> tuple[RegressionCase, ...]:
    root = tutorials_root_default()
    if not root.is_dir():
        return _KNOWN_CASES

    known_dirs = {c.case_dir for c in _KNOWN_CASES}
    discovered: list[RegressionCase] = list(_KNOWN_CASES)

    for ref_path in root.glob("**/regression/*.reference"):
        case_path = ref_path.parent.parent
        case_dir = str(case_path.relative_to(root))
        if case_dir not in known_dirs:
            discovered.append(
                RegressionCase(
                    case_dir=case_dir,
                    entry_name=None,
                    dicts=(),
                    reference_file=str(ref_path.relative_to(case_path)),
                )
            )
            known_dirs.add(case_dir)
            
    return tuple(discovered)


REGRESSION_CASES: tuple[RegressionCase, ...] = _discover_cases()
