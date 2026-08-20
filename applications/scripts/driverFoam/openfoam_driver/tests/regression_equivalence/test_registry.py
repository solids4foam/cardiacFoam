"""Registry guardrails for driverFOAM-owned regression-equivalence cases."""
from __future__ import annotations

import re
import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.specs.common import tutorials_root_default
from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES


def _alltest_entries() -> list[str]:
    script = tutorials_root_default() / "Alltest-regression"
    text = script.read_text()
    block = re.search(r"REGRESSION_TESTS=\((.*?)\)", text, re.DOTALL)
    assert block, "REGRESSION_TESTS array not found in Alltest-regression"
    return re.findall(r'"([^"]+)"', block.group(1))


def test_registry_cases_are_a_curated_subset_of_alltest_regression():
    registry_paths = {case.case_dir for case in REGRESSION_CASES}
    assert registry_paths <= set(_alltest_entries())


def test_every_case_dir_and_reference_exist_on_disk():
    root = tutorials_root_default()
    for case in REGRESSION_CASES:
        assert (root / case.case_dir).is_dir(), case.case_dir
        assert (root / case.case_dir / case.reference_file).is_file(), (
            f"{case.case_dir}/{case.reference_file}"
        )
        # assert (root / case.case_dir / case.regression_script).is_file(), (
        #     f"{case.case_dir}/{case.regression_script}"
        # )


def test_mapped_cases_have_entry_unmapped_have_none():
    for case in REGRESSION_CASES:
        if case.mapped:
            assert case.entry_name, case.case_dir
            assert case.drivers == ("strict", "generic"), case.case_dir
        else:
            assert case.entry_name is None, case.case_dir
            assert case.drivers == ("generic",), case.case_dir
