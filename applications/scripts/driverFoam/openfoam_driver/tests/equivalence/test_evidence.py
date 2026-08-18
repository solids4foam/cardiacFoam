"""L1 evidence capture: input-tree hashing and environment manifest."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


class TestInputTreeHash(unittest.TestCase):
    def test_hashes_every_file_by_relative_path(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import hash_input_tree
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "system").mkdir()
            (root / "system" / "controlDict").write_text("hello")
            hashes = hash_input_tree(root)
            self.assertEqual(
                hashes["system/controlDict"],
                "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824",
            )

    def test_excludes_volatile_paths(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import hash_input_tree
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "postProcessing").mkdir()
            (root / "postProcessing" / "out.dat").write_text("x")
            (root / "constant").mkdir()
            (root / "constant" / "keep").write_text("x")
            hashes = hash_input_tree(root)
            self.assertIn("constant/keep", hashes)
            self.assertNotIn("postProcessing/out.dat", hashes)

    def test_identical_trees_share_a_digest(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import (
            hash_input_tree, tree_digest,
        )
        digests = []
        for _ in range(2):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / "system").mkdir()
                (root / "system" / "controlDict").write_text("same")
                digests.append(tree_digest(hash_input_tree(root)))
        self.assertEqual(digests[0], digests[1])

    def test_differing_trees_differ_in_digest(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import (
            hash_input_tree, tree_digest,
        )
        digests = []
        for content in ("a", "b"):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / "system").mkdir()
                (root / "system" / "controlDict").write_text(content)
                digests.append(tree_digest(hash_input_tree(root)))
        self.assertNotEqual(digests[0], digests[1])


class TestEnvironmentManifest(unittest.TestCase):
    def test_records_openfoam_identity(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import environment_manifest
        manifest = environment_manifest({
            "WM_PROJECT_DIR": "/opt/openfoam",
            "WM_PROJECT_VERSION": "v2512",
            "WM_OPTIONS": "linux64GccDPInt32Opt",
            "IRRELEVANT": "ignored",
        })
        self.assertEqual(manifest["WM_PROJECT_DIR"], "/opt/openfoam")
        self.assertEqual(manifest["WM_PROJECT_VERSION"], "v2512")
        self.assertNotIn("IRRELEVANT", manifest)

    def test_absent_variables_are_none_not_missing(self) -> None:
        from openfoam_driver.tests.equivalence.evidence import environment_manifest
        manifest = environment_manifest({})
        self.assertIn("WM_PROJECT_DIR", manifest)
        self.assertIsNone(manifest["WM_PROJECT_DIR"])


if __name__ == "__main__":
    unittest.main()
