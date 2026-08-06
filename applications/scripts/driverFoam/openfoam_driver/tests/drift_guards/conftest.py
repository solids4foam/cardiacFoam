"""Module-level skip: all drift guards require the full cardiacFoam monorepo
tree (src/) to compare driverFOAM's Python representation against the C++ source of truth.

In a standalone driverFOAM clone the src/ tree is absent; pytest will skip
the entire test_*.py collection in this package.
"""
import pytest

from openfoam_driver.tests.conftest import monorepo_root

pytestmark = pytest.mark.skipif(
    monorepo_root is None,
    reason=(
        "Requires the full cardiacFoam monorepo tree (src/). "
        "Clone the full repository to enable drift guard tests."
    ),
)
