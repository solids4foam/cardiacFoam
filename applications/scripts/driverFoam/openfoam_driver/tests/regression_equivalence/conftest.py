"""Module-level skip: all regression-equivalence tests require the full
cardiacFoam monorepo tree (tutorials/ + applications/) to find the real case
directories, Alltest-regression script, and reference files.

In a standalone driverFOAM clone these files are absent; pytest will skip
the entire test_*.py collection in this package.
"""
import pytest

from openfoam_driver.tests.conftest import monorepo_root

pytestmark = pytest.mark.skipif(
    monorepo_root is None,
    reason=(
        "Requires the full cardiacFoam monorepo tree (tutorials/ + applications/). "
        "Clone the full repository to enable regression-equivalence tests."
    ),
)
