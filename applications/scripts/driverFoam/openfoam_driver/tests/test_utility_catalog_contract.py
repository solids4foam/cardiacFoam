"""Contract test: every utility directory has a manifest and the catalog is self-consistent.

When ``test_every_utility_dir_has_a_manifest`` fails, create a
``utility.manifest.toml`` in the missing utility directory.

When any other test fails, either the manifest TOML or the catalog loader has
drifted; fix the manifest file or the category constant in ``utility_catalog.py``.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.utility_catalog import (
    ALLOWED_CATEGORIES,
    MANIFEST_FILENAME,
    UTILITY_CATALOG,
    load_utility_manifests,
)

REPO_ROOT = Path(__file__).resolve().parents[5]
UTILITIES_ROOT = REPO_ROOT / "applications" / "utilities"


def _utility_dirs() -> list[Path]:
    """Return sorted list of direct subdirectories under UTILITIES_ROOT."""
    return sorted(p for p in UTILITIES_ROOT.iterdir() if p.is_dir())


class TestUtilityCatalogContract(unittest.TestCase):
    def test_every_utility_dir_has_a_manifest(self) -> None:
        """Every directory under applications/utilities/ must have a manifest file."""
        missing = [
            d.name
            for d in _utility_dirs()
            if not (d / MANIFEST_FILENAME).exists()
        ]
        self.assertEqual(
            missing,
            [],
            f"Utility directories missing {MANIFEST_FILENAME}: {missing}",
        )

    def test_manifest_name_matches_directory(self) -> None:
        """Each loaded manifest's 'name' field must equal its parent directory name."""
        mismatches = []
        for name, manifest in UTILITY_CATALOG.items():
            dir_name = manifest.source_path.parent.name
            if manifest.name != dir_name:
                mismatches.append(
                    f"{manifest.source_path}: name={manifest.name!r} "
                    f"vs dir={dir_name!r}"
                )
        self.assertEqual(
            mismatches,
            [],
            f"name/directory mismatches: {mismatches}",
        )

    def test_no_duplicate_manifest_names(self) -> None:
        """The catalog size must equal the number of utility directories that have manifests."""
        dirs_with_manifests = [
            d for d in _utility_dirs() if (d / MANIFEST_FILENAME).exists()
        ]
        self.assertEqual(
            len(UTILITY_CATALOG),
            len(dirs_with_manifests),
            f"Catalog has {len(UTILITY_CATALOG)} entries but "
            f"{len(dirs_with_manifests)} manifest files exist; "
            "check for duplicate 'name' values across manifests.",
        )

    def test_required_fields_present(self) -> None:
        """Every manifest must have non-empty name, description, and category."""
        violations = []
        for name, manifest in UTILITY_CATALOG.items():
            if not manifest.name:
                violations.append(f"{name}: 'name' is empty")
            if not manifest.description:
                violations.append(f"{name}: 'description' is empty")
            if not manifest.category:
                violations.append(f"{name}: 'category' is empty")
        self.assertEqual(
            violations,
            [],
            f"Required-field violations: {violations}",
        )

    def test_categories_in_allowed_set(self) -> None:
        """Every manifest's category must be one of ALLOWED_CATEGORIES."""
        bad = [
            f"{name}: {manifest.category!r}"
            for name, manifest in UTILITY_CATALOG.items()
            if manifest.category not in ALLOWED_CATEGORIES
        ]
        self.assertEqual(
            bad,
            [],
            f"Categories not in ALLOWED_CATEGORIES {sorted(ALLOWED_CATEGORIES)}: {bad}",
        )


if __name__ == "__main__":
    unittest.main()
