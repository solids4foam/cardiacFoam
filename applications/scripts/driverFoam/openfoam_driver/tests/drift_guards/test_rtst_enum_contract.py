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
#     test_rtst_enum_contract
#
# Description
#     Tests rtst enum contract logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Contract test: catalogue `enum_values` match registered runtime-selection
types in `src/`.

Drift is reported in two directions:

* a catalogue enum's values diverge from the RTST it claims to expose, or
* an RTST exists in code with no driver-facing catalogue mapping and no
  entry in `INTERNAL_RTST_ALLOWLIST` (i.e. an undocumented hierarchy).

The mapping `RTST_BY_DRIVER_PATH` below names the RTST base each catalogue
enum corresponds to. New bases or new catalogue enums must be wired up
here explicitly — the test refuses to silently pass on unmapped entries.
"""

from __future__ import annotations

import unittest
from pathlib import Path
import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo

from openfoam_driver.scripts._rtst_scanner import (
    iter_catalogue_enums,
    scan_rtst_registrations,
)

REPO_ROOT = Path(__file__).resolve().parents[6]
SRC_ROOT = REPO_ROOT / "src"


# Catalogue driver_path → (RTST base class, comparison mode).
#
# Modes:
#   "strict" — catalogue enum_values must equal the set of registered names.
#   "subset" — catalogue enum_values must be a subset of registered names
#              (used where one driver field deliberately exposes a curated
#              subset of an RTST, e.g. Purkinje only uses a subset of the
#              ionicModel hierarchy).
RTST_BY_DRIVER_PATH: dict[str, tuple[str, str]] = {
    "type": ("physicsModel", "strict"),
    "$ELECTRO_MODEL_COEFFS.ionicModel": ("ionicModel", "strict"),
    "$ELECTRO_MODEL_COEFFS.electrophysicsAdvanceScheme": (
        "electrophysicsAdvanceScheme",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.verificationModel.type": (
        "electroVerificationModel",
        "polymorphic",
    ),
    # The Purkinje graph and the PVJ coupling each have their own verifier
    # table, separate from the myocardium-side electroVerificationModel.
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".purkinjeGraphModelCoeffs.verificationModel.type": (
        "graphVerificationModel",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.verificationModel.type": (
        "couplingVerificationModel",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgSolver": (
        "ecgSolver",
        "strict",
    ),
    # ecgVerificationModel.H:73 declares a runtime selection table; the enum
    # values are registered verifier type names. Polymorphic because the
    # bath/pseudo/eikonal verifiers register against the same base.
    "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.ecgVerificationModel": (
        "ecgVerificationModel",
        "polymorphic",
    ),
    "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.coupling.electroDomainCoupler": (
        "electroDomainCoupler",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".purkinjeGraphModelCoeffs.conductionSystemSolver": (
        "conductionSystemSolver",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".purkinjeGraphModelCoeffs.ionicModel": ("ionicModel", "subset"),
    "$ELECTRO_MODEL_COEFFS.activeTensionModel": (
        "activeTensionModel",
        "strict",
    ),
    "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.electroDomainCoupler": (
        "electroDomainCoupler",
        "strict",
    ),
}


# Catalogue enums that are NOT runtime-selection tables in the C++ code
# (their values come from manual enums or external libraries). Listing
# them here suppresses the "unmapped catalogue enum" warning.
NON_RTST_DRIVER_PATHS: frozenset[str] = frozenset({
    "$ELECTRO_MODEL_COEFFS.solutionAlgorithm",       # implicit/explicit
    "$ELECTRO_MODEL_COEFFS.tissue",                  # tissue enum
    "$ELECTRO_MODEL_COEFFS.sex",                     # biological-sex selector
    "$ELECTRO_MODEL_COEFFS.dimension",               # 1D/2D/3D
    # Same 1D/2D/3D selector at Purkinje scope: ionicSelector maps the word
    # to a flag, there is no runtime table.
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".purkinjeGraphModelCoeffs.dimension",
    "$ELECTRO_MODEL_COEFFS.solver",                  # OpenFOAM ODESolver
    "$ELECTRO_MODEL_COEFFS.couplingSignal",
    # timeCouplingScheme is a plain lookupOrDefault<word> compared against two
    # literals in myocardiumDomain.C:359-365, with a FatalError otherwise --
    # no runtime table, no registered classes.
    "$ELECTRO_MODEL_COEFFS.timeCouplingScheme",
    "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.couplingMode",
    "$ELECTRO_MODEL_COEFFS.ecgDomains.<name>.manufactured.dimension",
    "$ELECTRO_MODEL_COEFFS.verificationModel.fdaBathVariant",
    "$ELECTRO_MODEL_COEFFS.manufacturedBidomain.fdaBathVariant",
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".purkinjeGraphModelCoeffs.tissue",
    # pvjKernel selects the PVJ spatial-spreading kernel (uniform / ...) via a
    # plain lookupOrDefault<word> in pvjCoupler.C, not a runtime table.
    "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjKernel",
    # pvjCouplingScheme selects explicit vs split implicit source assembly by
    # a plain word in reactionDiffusionPvjCoupler.C, not a runtime table.
    "$ELECTRO_MODEL_COEFFS.domainCouplings.<name>.pvjCouplingScheme",
    # The user-facing `myocardiumSolver` value is a UNION of the
    # `electroModel` RTST and the `myocardiumSolver` RTST, plus the
    # eikonalMyocardiumDomain factory — no single base to check against.
    "myocardiumSolver",
    # The `conductionSystemDomain` selector has a single accepted value
    # (`purkinjeGraphModel`); the value is not registered in any RTST.
    "$ELECTRO_MODEL_COEFFS.conductionNetworkDomains.<name>"
    ".conductionSystemDomain",
    # batchedIntegrator selects the ODE integrator for compactBatched GPU
    # models (euler / rushLarsen / rushLarsenHeun). Not an RTST — it is an
    # internal solver-algorithm enum with no separate runtime table.
    "$ELECTRO_MODEL_COEFFS.batchedIntegrator",
    # ionicHeterogeneity enums are plain dictionary words parsed directly by
    # ionicHeterogeneity.C / ionicModel.C (lookupOrDefault), not RTST types.
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mode",            # transmuralBands
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode",  # blend/hard
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing",       # smoothstep
    # extracellular bath-assembly enums parsed directly by
    # extracellularPotentialDomain.C (lookupOrDefault<word>), not RTST types.
    "$ELECTRO_MODEL_COEFFS.bathPotentialDomain.interfaceConductivityInterpolation",
    # conductivitySource selects how the myocardial conductivity tensor is
    # obtained (uniform dict literal / field read from disk). It is parsed as
    # a plain word by conductivityFieldIO.C (lookup + manual branch), not an
    # RTST type.
    "$ELECTRO_MODEL_COEFFS.conductivitySource",
    # fdaBathVariant (groundElectrode/electrodePair) selects the outer bath
    # boundary-condition layout for the FDA bath-bidomain manufactured
    # verifier. It is a plain NamedEnum word read by
    # manufacturedFDABathBidomainVerifier.C, not an RTST runtime-selection
    # table.
    "$ELECTRO_MODEL_COEFFS.verificationModel.fdaBathVariant",
})


# RTST bases that legitimately have no driver-facing catalogue enum.
# Either they are consumed indirectly (the user does not type the value
# directly into a dictionary) or they have only one registered type and
# selection is implicit. Listing a base here documents the omission.
INTERNAL_RTST_ALLOWLIST: frozenset[str] = frozenset({
    # `electroModel` and `myocardiumSolver` are the lower-level tables the
    # `myocardiumSolver` dict key resolves through; values appear in the
    # union enum, not in their own dedicated catalogue entry.
    "electroModel",
    "myocardiumSolver",
    # `electroMechanicalModel` is selected via the top-level
    # PHYSICS_PROPERTY_ENTRIES `type` enum (the single registered
    # `sequentialElectroMechanical` is the only option once that path is
    # taken; no driver-side enum needed yet).
    "electroMechanicalModel",
    # `ecgVerificationModel` registers verifier types but the catalogue
    # currently exposes only `manufactured.enabled` (boolean) — no
    # ecg-side `type` selector yet.
    "ecgVerificationModel",
    # Only one bath-potential implementation exists today. Users select it
    # by adding bathPotentialDomain, not by typing an electroStateDomain type.
    "electroStateDomain",
    # Like ecgVerificationModel, there is no catalogue type selector for this yet.
    "electromechanicalVerificationModel",
    "eikonalVerificationModel",
    # Coupling- and graph-side manufactured verifiers (1D-3D Purkinje work).
    # Like the other *VerificationModel hierarchies, selected via the
    # verificationModel subdict, not a dedicated catalogue type enum.
    "couplingVerificationModel",
    "graphVerificationModel",
})


@skip_without_monorepo
class TestRtstEnumContract(unittest.TestCase):
    """All bases discovered by the scanner must be addressed by either
    a `RTST_BY_DRIVER_PATH` mapping or an `INTERNAL_RTST_ALLOWLIST`
    entry. The catalogue's enum_values must agree with the registered
    set for every mapped path.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.regs = scan_rtst_registrations(SRC_ROOT)
        cls.catalogue_enums = list(iter_catalogue_enums())
        cls.mapped_bases = {
            base for base, _ in RTST_BY_DRIVER_PATH.values()
        }

    def test_every_registered_base_is_mapped_or_allowlisted(self) -> None:
        """No silent drift: a base that appears in `addToRunTimeSelectionTable`
        must either drive a catalogue enum or be explicitly skipped.
        """
        unmapped = []
        for base in self.regs:
            if base in self.mapped_bases:
                continue
            if base in INTERNAL_RTST_ALLOWLIST:
                continue
            unmapped.append(base)
        self.assertEqual(
            unmapped,
            [],
            "RTST bases registered in src/ but neither mapped to a "
            f"catalogue enum nor allowlisted: {sorted(unmapped)}. "
            "Add a RTST_BY_DRIVER_PATH or INTERNAL_RTST_ALLOWLIST entry.",
        )

    def test_every_catalogue_enum_is_classified(self) -> None:
        """Every catalogue enum must either map to an RTST base or be
        listed as a non-RTST enum. Forces deliberate handling of new
        enum entries.
        """
        unclassified = []
        for c in self.catalogue_enums:
            if c.driver_path in RTST_BY_DRIVER_PATH:
                continue
            if c.driver_path in NON_RTST_DRIVER_PATHS:
                continue
            unclassified.append(c.driver_path)
        self.assertEqual(
            unclassified,
            [],
            "Catalogue enums with no RTST_BY_DRIVER_PATH or "
            f"NON_RTST_DRIVER_PATHS entry: {unclassified}",
        )

    def test_mapped_catalogue_enums_agree_with_registrations(self) -> None:
        for c in self.catalogue_enums:
            if c.driver_path not in RTST_BY_DRIVER_PATH:
                continue
            base, mode = RTST_BY_DRIVER_PATH[c.driver_path]
            with self.subTest(driver_path=c.driver_path, base=base, mode=mode):
                self.assertIn(
                    base,
                    self.regs,
                    f"mapping rule references base {base!r} but no "
                    f"addToRunTimeSelectionTable({base}, ...) was found",
                )
                registered = set(self.regs[base])
                catalogue_set = set(c.enum_values)
                if mode == "strict":
                    self.assertEqual(
                        catalogue_set,
                        registered,
                        f"{c.driver_path}: catalogue vs RTST {base} "
                        f"differ. "
                        f"in catalogue only: {sorted(catalogue_set - registered)}; "
                        f"in RTST only: {sorted(registered - catalogue_set)}.",
                    )
                elif mode == "subset":
                    extra = catalogue_set - registered
                    self.assertEqual(
                        extra,
                        set(),
                        f"{c.driver_path} (subset mode): catalogue lists "
                        f"values not registered in {base}: {sorted(extra)}",
                    )
                elif mode == "polymorphic":
                    missing_from_schema = registered - catalogue_set
                    self.assertEqual(
                        missing_from_schema,
                        set(),
                        f"{c.driver_path} (polymorphic mode): RTST {base} has "
                        f"items not in python schema catalogue: {sorted(missing_from_schema)}",
                    )
                else:
                    self.fail(f"unknown comparison mode: {mode!r}")


if __name__ == "__main__":
    unittest.main()
