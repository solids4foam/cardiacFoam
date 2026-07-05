"""Agent dict-layer round-trip idempotence for mapped regressions (solver-free)."""
from __future__ import annotations

import difflib

import pytest

from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES
from openfoam_driver.tests.regression_equivalence.round_trip import (
    electro_build_parse_fixpoint,
)

MAPPED = [c for c in REGRESSION_CASES if c.mapped]


@pytest.mark.parametrize("case", MAPPED, ids=lambda c: c.entry_name)
def test_electro_build_parse_is_idempotent(case):
    once, twice = electro_build_parse_fixpoint(case)
    if once != twice:
        diff = "".join(
            difflib.unified_diff(
                once.splitlines(keepends=True),
                twice.splitlines(keepends=True),
                fromfile="build∘parse", tofile="build∘parse∘build∘parse",
            )
        )
        pytest.fail(f"{case.entry_name}: build∘parse not idempotent\n{diff}")
