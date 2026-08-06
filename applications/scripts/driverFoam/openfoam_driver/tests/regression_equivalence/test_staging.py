"""Agent addressability of every regression case (solver-free)."""
from __future__ import annotations

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES
from openfoam_driver.tests.regression_equivalence.staging import (
    resolve_generic,
    resolve_strict,
)

ADDRESSABLE = [c for c in REGRESSION_CASES if c.generic_addressable]
NOT_ADDRESSABLE = [c for c in REGRESSION_CASES if not c.generic_addressable]
MAPPED = [c for c in REGRESSION_CASES if c.mapped]


@pytest.mark.parametrize("case", ADDRESSABLE, ids=lambda c: c.case_dir)
def test_generic_path_resolves(case):
    resolution = resolve_generic(case)
    assert resolution["resolution"] == "case_folder", case.case_dir
    assert resolution["is_runnable"] is True, case.case_dir


@pytest.mark.parametrize("case", NOT_ADDRESSABLE, ids=lambda c: c.case_dir)
def test_non_addressable_case_is_not_discoverable(case):
    # Documents a real limitation: the agent's case discovery does not recognize
    # this layout, so it cannot be driven generically.
    with pytest.raises(KeyError, match="Unknown entry"):
        resolve_generic(case)


@pytest.mark.parametrize("case", MAPPED, ids=lambda c: c.entry_name)
def test_mapped_entry_resolves_registered(case):
    resolution = resolve_strict(case)
    assert resolution["resolution"] == "registered", case.entry_name
    assert resolution["is_runnable"] is True, case.entry_name
