"""manufacturedBathBidomain must apply cleanly to its own checked-in case.

The bath-bidomain tutorial switches between two mutually exclusive boundary
variants. ``groundPatches`` and ``surfaceCurrentPatches`` may not both claim
the same patch (``extracellularPotentialDomain.C`` rejects it), so switching
variants must *remove* the other variant's patch entry, not merely add its own.

The committed ``electroProperties`` records whichever variant ran last -- the
shared ``case_root`` that entry-based sweeps mutate in place. So applying the
tutorial always has to reconcile against a case that may be in either state,
and the removal half of that reconciliation has to actually work.

It did not. Both cleanups reached for ``remove_electro_property_dict``, a
*block* remover, to delete ``xMin`` -- which is a scalar entry (``xMin -0.01;``),
not a block. The remover matched the line, scanned forward for an opening
brace, found none, and raised ``KeyError: "Dictionary 'xMin' has no opening
brace"``. ``missing_ok=True`` did not help: it wraps only scope resolution, so
a key that is present but the wrong *shape* passes straight through.

Nothing about the case was invalid. ``groundElectrode`` is a valid C++ variant
and the dictionary is valid OpenFOAM; the fault was entirely in driverFOAM's
mutation layer, raised before any solver ran.

No existing test caught it: the tet suite passes ``mesh_family="tet"`` and
builds against ``tmp_path``, never exercising the default hex path against the
real tracked case.
"""

from __future__ import annotations

import re
import shutil
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.runtime.models import CaseConfig
from openfoam_driver.plugins.cardiacfoam.tutorials import manufactured_bath_bidomain as tutorial
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


def _patch_entries(text: str, block: str) -> dict[str, str]:
    """Return the scalar entries inside one named patch sub-dictionary."""
    match = re.search(rf"{block}\s*\{{(.*?)\}}", text, re.S)
    assert match is not None, f"{block} block not found"
    return dict(re.findall(r"(\w+)\s+(\S+?);", match.group(1)))


@skip_without_monorepo
class TestBathBidomainVariantApply(unittest.TestCase):
    """Both variants must apply to a pristine copy of the tracked case."""

    def setUp(self) -> None:
        source = (
            monorepo_root / "tutorials" / "manufacturedSolutions" / "bathBidomain"
        )
        if not source.is_dir():
            self.skipTest(f"tracked bathBidomain case not found at {source}")
        self.tutorials_root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        destination = self.tutorials_root / "manufacturedSolutions" / "bathBidomain"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, destination)
        self.case_root = destination
        self.electro = destination / "constant" / "electroProperties"

    def _apply(self, variant: str) -> str:
        spec = tutorial.make_spec(
            tutorials_root=self.tutorials_root, fda_bath_variant=variant
        )
        spec.apply_case(spec.case_root, spec.build_cases()[0])
        return self.electro.read_text()

    def test_default_variant_applies_to_the_tracked_case(self) -> None:
        """The regression: make_spec() defaults, pristine tracked case."""
        spec = tutorial.make_spec(tutorials_root=self.tutorials_root)
        spec.apply_case(spec.case_root, spec.build_cases()[0])

    def test_ground_electrode_removes_the_electrode_pair_patch(self) -> None:
        text = self._apply("groundElectrode")

        ground = _patch_entries(text, "groundPatches")
        surface = _patch_entries(text, "surfaceCurrentPatches")

        self.assertIn("xMin", ground, "groundElectrode must claim xMin")
        self.assertNotIn(
            "xMin",
            surface,
            "xMin may not be claimed by both patch lists -- "
            "extracellularPotentialDomain.C rejects that",
        )

    def test_electrode_pair_removes_the_ground_patch(self) -> None:
        # Force the opposite starting state first, so the cleanup has real work.
        self._apply("groundElectrode")
        text = self._apply("electrodePair")

        ground = _patch_entries(text, "groundPatches")
        surface = _patch_entries(text, "surfaceCurrentPatches")

        self.assertIn("xMin", surface, "electrodePair must claim xMin")
        self.assertNotIn("xMin", ground, "xMin may not be claimed by both lists")

    def test_variant_switching_is_idempotent_and_reversible(self) -> None:
        """Applying a variant must not depend on which variant preceded it.

        This is the mutation-totality invariant in miniature: the shared
        case_root means a key set on one path must be removed on the other, or
        the result depends on invocation history rather than on the request.
        """
        from_clean = self._apply("groundElectrode")
        after_round_trip = (self._apply("electrodePair"), self._apply("groundElectrode"))[-1]

        self.assertEqual(
            _patch_entries(from_clean, "groundPatches"),
            _patch_entries(after_round_trip, "groundPatches"),
        )
        self.assertEqual(
            _patch_entries(from_clean, "surfaceCurrentPatches"),
            _patch_entries(after_round_trip, "surfaceCurrentPatches"),
        )


if __name__ == "__main__":
    unittest.main()
