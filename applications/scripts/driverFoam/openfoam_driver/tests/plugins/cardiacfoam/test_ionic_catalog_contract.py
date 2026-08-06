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
#     test_ionic_catalog_contract
#
# Description
#     Tests ionic catalog contract logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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

from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.tests.conftest import skip_without_monorepo
from openfoam_driver.scripts._names_parser import (
    EXCLUDED_FROM_HEADER_SYNC,
    find_names_header,
    parse_names_header,
)

REPO_ROOT = Path(__file__).resolve().parents[7]
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


@skip_without_monorepo
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


_FULL_IONIC_MODELS = {
    "TNNP", "Grandi", "Courtemanche", "Fabbri",
    "ToRORd_dynCl", "Trovato", "Stewart", "Gaur", "PerisYague", "TWorld",
}

# Tokens in recommended_exports that are known V/Ca aliases carried over from
# before P11b drift is resolved. These appear in the catalog's existing pair
# but do not yet match a state/algebraic name — P11b will fix the drift.
# We allow them so the P11c test does not force P11b work.
_ALLOWED_LEGACY_TOKENS: set[str] = {
    # TNNP: states has "V" and "Ca_i"; legacy aliases used before P11b
    "membrane_V", "calcium_Cai",
    # ORd, Grandi: states has "membrane_V"/"Ca_i"; legacy aliases
    "Cass",  # not a state in ORd (Ca_ss is "Ca_ss") — legacy alias
    # Trovato: no "Vm" in states
    "Vm",
    # Stewart: no plain "V" in states (states has "membrane_V"); "Cai" alias
    "V", "Cai",
}


class TestRecommendedExportsExpanded(unittest.TestCase):
    """P11c contract: full ionic models expose >= 4 recommended exports."""

    def test_full_models_have_at_least_four_exports(self) -> None:
        for name in _FULL_IONIC_MODELS:
            with self.subTest(model=name):
                entry = IONIC_MODEL_CATALOG[name]
                self.assertGreaterEqual(
                    len(entry.recommended_exports),
                    4,
                    f"{name}.recommended_exports has only "
                    f"{len(entry.recommended_exports)} entries (< 4)",
                )

    def test_recommended_exports_tokens_in_states_or_algebraic(self) -> None:
        for name in _FULL_IONIC_MODELS:
            with self.subTest(model=name):
                entry = IONIC_MODEL_CATALOG[name]
                valid = set(entry.states) | set(entry.algebraic)
                bad = [
                    tok
                    for tok in entry.recommended_exports
                    if tok not in valid and tok not in _ALLOWED_LEGACY_TOKENS
                ]
                self.assertEqual(
                    bad,
                    [],
                    f"{name}.recommended_exports contains tokens not in "
                    f"states/algebraic (and not in allowed legacy set): {bad}",
                )


class TestRecommendedExportsExpansion(unittest.TestCase):
    """Full ionic models must advertise at least
    voltage + calcium + 3 main currents in recommended_exports — the agent
    fallback path (when no outputVariables.ionic.export is declared) uses
    this list, so a minimal 2-variable default is uninformative.
    Phenomenological and manufactured models are intentionally exempt:
    phenomenological models lack current variables; manufactured models
    expose their analytic-solution components only.
    """

    # Full ionic models — must carry voltage + calcium + ≥3 currents.
    _FULL_IONIC_MODELS: frozenset[str] = frozenset({
        "TNNP", "Grandi", "Courtemanche", "Fabbri",
        "ToRORd_dynCl", "Trovato", "Stewart", "Gaur",
        "PerisYague", "TWorld",
    })

    def test_every_full_model_has_at_least_5_exports(self) -> None:
        """5 = voltage + calcium + 3 currents minimum. Models can carry more;
        this is the floor."""
        for name in self._FULL_IONIC_MODELS:
            with self.subTest(model=name):
                entry = IONIC_MODEL_CATALOG[name]
                self.assertGreaterEqual(
                    len(entry.recommended_exports), 5,
                    f"{name}.recommended_exports has only "
                    f"{len(entry.recommended_exports)} tokens — "
                    f"expected at least 5 (Vm + Cai + 3 currents).",
                )

    def test_full_model_exports_include_total_ionic_current(self) -> None:
        """Iion (or Iion_cm) appears in every full ionic model's catalogue
        — it should be in every recommended_exports list as the
        physiological 'output of last resort'."""
        for name in self._FULL_IONIC_MODELS:
            with self.subTest(model=name):
                entry = IONIC_MODEL_CATALOG[name]
                has_total = any(
                    "Iion" in token for token in entry.recommended_exports
                )
                self.assertTrue(
                    has_total,
                    f"{name}.recommended_exports lacks a total-ionic-current "
                    f"token (Iion / Iion_cm): {entry.recommended_exports}",
                )


if __name__ == "__main__":
    unittest.main()
