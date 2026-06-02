"""Audit test: ionic_model_catalog.py Role-A fields vs C++ Names.H enums.

This test is the regression guard described in plan §11.2 (2026-05-19-driverfoam-
agentic-integration.md).  It directly parses each model's ``*_Names.H`` with an
inline regex, extracts the three enum bodies in declaration order, and asserts
exact equality against the catalogue's ``states``, ``algebraic``, and
``constants`` tuples.

Design choices
--------------
* Manufactured-family models (monodomainFDAManufactured, bidomainFDAManufactured,
  bathBidomainFDAManufactured) are excluded — their catalogue entries carry
  user-facing semantic labels, not raw C++ enum tokens.  See §3d-3 in the plan.
* Only Role-A fields (``states``, ``algebraic``, ``constants``) are checked.
  Role-B metadata (``compatible_solvers``, ``species``, ``description``, etc.)
  is hand-curated and not derivable from C++.
* The enum-name matching is intentionally liberal: we look for any enum whose
  name contains "STATE", "ALGEBRAIC", or "CONSTANT" (case-insensitive) to handle
  variants like ``STATES_INDEX``, ``STATE_INDEX``, ``CONSTANT_INDEX``,
  ``CONSTANTS_INDEX``.
* ``NUM_*`` sentinels are stripped before comparison.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG
from openfoam_driver.scripts._names_parser import EXCLUDED_FROM_HEADER_SYNC

# ---------------------------------------------------------------------------
# Repository layout
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[5]
IONIC_MODELS_DIR = REPO_ROOT / "src" / "ionicModels"

# ---------------------------------------------------------------------------
# Models to audit
# ---------------------------------------------------------------------------
_MANUFACTURED_FAMILY: frozenset[str] = frozenset(
    {
        "monodomainFDAManufactured",
        "bidomainFDAManufactured",
        "bathBidomainFDAManufactured",
    }
)

NON_MANUFACTURED_MODELS: list[str] = [
    name
    for name in IONIC_MODEL_CATALOG
    if name not in _MANUFACTURED_FAMILY
]

# Models audited against C++ Names.H headers: excludes manufactured family
# and any models explicitly excluded from header sync (e.g. compactBatched
# GPU models that reuse the parent CPU model's header).
HEADER_AUDITED_MODELS: list[str] = [
    name for name in NON_MANUFACTURED_MODELS
    if name not in EXCLUDED_FROM_HEADER_SYNC
]

# ---------------------------------------------------------------------------
# Inline parser (no dependency on the regenerate-script infrastructure)
# ---------------------------------------------------------------------------

_BLOCK_COMMENT_RE = re.compile(r"/\*.*?\*/", re.DOTALL)
_LINE_COMMENT_RE = re.compile(r"//[^\n]*")
_ENUM_HEADER_RE = re.compile(r"enum\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", re.MULTILINE)


def _strip_comments(text: str) -> str:
    text = _BLOCK_COMMENT_RE.sub("", text)
    text = _LINE_COMMENT_RE.sub("", text)
    return text


def _is_sentinel(token: str) -> bool:
    upper = token.upper()
    return upper.startswith("NUM_") or "_NUM_" in upper


def _parse_enum_body(body: str) -> tuple[str, ...]:
    body = _strip_comments(body)
    items: list[str] = []
    for raw in body.split(","):
        token = raw.strip()
        if not token:
            continue
        if "=" in token:
            token = token.split("=", 1)[0].strip()
        if not token:
            continue
        if _is_sentinel(token):
            continue
        items.append(token)
    return tuple(items)


def _extract_enum(text: str, keyword: str) -> tuple[str, ...] | None:
    """Return tokens from the first enum whose name contains *keyword*."""
    clean = _strip_comments(text)
    kw_up = keyword.upper()
    for m in _ENUM_HEADER_RE.finditer(clean):
        if kw_up not in m.group(1).upper():
            continue
        start = m.end()
        depth = 1
        i = start
        while i < len(clean) and depth > 0:
            ch = clean[i]
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
            i += 1
        if depth != 0:
            return None
        return _parse_enum_body(clean[start : i - 1])
    return None


def _find_names_header(model_dir: Path) -> Path | None:
    if not model_dir.is_dir():
        return None
    candidates = [p for p in model_dir.glob("*Names.H") if "lnInclude" not in p.parts]
    return candidates[0] if candidates else None


def _parse_model(model: str) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    header = _find_names_header(IONIC_MODELS_DIR / model)
    if header is None:
        raise FileNotFoundError(f"No *Names.H found under {IONIC_MODELS_DIR / model}")
    text = header.read_text(encoding="utf-8", errors="replace")
    states = _extract_enum(text, "STATE") or ()
    algebraic = _extract_enum(text, "ALGEBRAIC") or ()
    constants = _extract_enum(text, "CONSTANT") or ()
    return states, algebraic, constants


# ---------------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------------


class TestIonicCatalogAudit(unittest.TestCase):
    """Parametrised audit: each non-manufactured model gets a subTest per field."""

    def _assert_field_matches(
        self,
        model: str,
        field: str,
        cpp_tokens: tuple[str, ...],
        cat_tokens: tuple[str, ...],
    ) -> None:
        cpp_set = set(cpp_tokens)
        cat_set = set(cat_tokens)
        extra_in_cpp = [t for t in cpp_tokens if t not in cat_set]
        extra_in_cat = [t for t in cat_tokens if t not in cpp_set]
        self.assertEqual(
            cpp_tokens,
            cat_tokens,
            f"{model}.{field} mismatch between C++ Names.H and catalog.\n"
            f"  Extra in C++ (missing from catalog): {extra_in_cpp}\n"
            f"  Extra in catalog (stale): {extra_in_cat}\n"
            f"  Ordering also compared — full diff visible above.",
        )

    def test_states_audit(self) -> None:
        for model in HEADER_AUDITED_MODELS:
            with self.subTest(model=model):
                cpp_states, _, _ = _parse_model(model)
                cat_states = IONIC_MODEL_CATALOG[model].states
                self._assert_field_matches(model, "states", cpp_states, cat_states)

    def test_algebraic_audit(self) -> None:
        for model in HEADER_AUDITED_MODELS:
            with self.subTest(model=model):
                _, cpp_algebraic, _ = _parse_model(model)
                cat_algebraic = IONIC_MODEL_CATALOG[model].algebraic
                self._assert_field_matches(model, "algebraic", cpp_algebraic, cat_algebraic)

    def test_constants_audit(self) -> None:
        for model in HEADER_AUDITED_MODELS:
            with self.subTest(model=model):
                _, _, cpp_constants = _parse_model(model)
                cat_constants = IONIC_MODEL_CATALOG[model].constants
                self._assert_field_matches(model, "constants", cpp_constants, cat_constants)

    def test_manufactured_family_not_audited(self) -> None:
        """Sanity-check: manufactured models are excluded from the audit list."""
        for name in _MANUFACTURED_FAMILY:
            self.assertNotIn(name, NON_MANUFACTURED_MODELS)

    def test_all_24_non_manufactured_models_present(self) -> None:
        """All 24 expected non-manufactured models appear in the audit list.

        12 CPU models + 12 compactBatched GPU models.
        ORd removed (C++ implementation deleted as non-functional).
        """
        expected = {
            # CPU models (12)
            "AlievPanfilov",
            "BuenoOrovio",
            "Courtemanche",
            "Fabbri",
            "Gaur",
            "Grandi",
            "PerisYague",
            "Stewart",
            "TNNP",
            "ToRORd_dynCl",
            "Trovato",
            "TWorld",
            # GPU compactBatched models (11)
            "AlievPanfilovcompactBatched",
            "BuenoOroviocompactBatched",
            "CourtemanchecompactBatched",
            "FabbricompactBatched",
            "GaurcompactBatched",
            "GrandicompactBatched",
            "PerisYaguecompactBatched",
            "StewartcompactBatched",
            "TNNPcompactBatched",
            "ToRORd_dynClcompactBatched",
            "TrovatocompactBatched",
            "TWorldcompactBatched",
        }
        self.assertEqual(
            expected,
            set(NON_MANUFACTURED_MODELS),
            "Catalog has unexpected models or is missing expected ones.",
        )


if __name__ == "__main__":
    unittest.main()
