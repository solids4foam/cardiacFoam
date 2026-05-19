"""Contract test: `ionic_model_catalog.py` matches each model's `*_Names.H`.

When this fails, the catalogue has drifted from the C++ source of truth. Fix
by running:

    python applications/scripts/driverFoam/scripts/regenerate-ionic-catalog.py

Models in `EXCLUDED_FROM_HEADER_SYNC` (the FDA manufactured family) are
intentionally documented with user-facing semantic labels rather than raw
enum identifiers and are not checked here.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.scripts._names_parser import (
    EXCLUDED_FROM_HEADER_SYNC,
    find_names_header,
    parse_names_header,
)

REPO_ROOT = Path(__file__).resolve().parents[5]
IONIC_MODELS_DIR = REPO_ROOT / "src" / "ionicModels"

_REGEN_HINT = (
    "run: python applications/scripts/driverFoam/scripts/"
    "regenerate-ionic-catalog.py"
)


def _models_to_check() -> list[str]:
    return [
        name
        for name in IONIC_MODEL_CATALOG
        if name not in EXCLUDED_FROM_HEADER_SYNC
    ]


class TestIonicCatalogContract(unittest.TestCase):
    def test_every_checked_model_has_a_names_header(self) -> None:
        missing = []
        for name in _models_to_check():
            header = find_names_header(IONIC_MODELS_DIR / name)
            if header is None:
                missing.append(name)
        self.assertEqual(
            missing,
            [],
            f"models with no *_Names.H header (and not in "
            f"EXCLUDED_FROM_HEADER_SYNC): {missing}",
        )

    def test_states_match_header(self) -> None:
        for name in _models_to_check():
            with self.subTest(model=name):
                header = find_names_header(IONIC_MODELS_DIR / name)
                self.assertIsNotNone(header, f"no header for {name}")
                assert header is not None  # for type-checker
                parsed = parse_names_header(header)
                entry = IONIC_MODEL_CATALOG[name]
                self.assertEqual(
                    entry.states,
                    parsed.states,
                    f"{name}.states drift from {header.name}; {_REGEN_HINT}",
                )

    def test_algebraic_match_header(self) -> None:
        for name in _models_to_check():
            with self.subTest(model=name):
                header = find_names_header(IONIC_MODELS_DIR / name)
                assert header is not None
                parsed = parse_names_header(header)
                entry = IONIC_MODEL_CATALOG[name]
                self.assertEqual(
                    entry.algebraic,
                    parsed.algebraic,
                    f"{name}.algebraic drift from {header.name}; {_REGEN_HINT}",
                )

    def test_constants_match_header(self) -> None:
        for name in _models_to_check():
            with self.subTest(model=name):
                header = find_names_header(IONIC_MODELS_DIR / name)
                assert header is not None
                parsed = parse_names_header(header)
                entry = IONIC_MODEL_CATALOG[name]
                self.assertEqual(
                    entry.constants,
                    parsed.constants,
                    f"{name}.constants drift from {header.name}; {_REGEN_HINT}",
                )


if __name__ == "__main__":
    unittest.main()
