"""Module-level skip: verification tests check that physical tutorial case
directories exist on disk and that sweep JSON files reference real paths.

In a standalone driverFOAM clone the tutorials/ tree is absent; pytest will
skip the entire test_*.py collection in this package.
"""
import pytest

from openfoam_driver.tests.conftest import monorepo_root

pytestmark = pytest.mark.skipif(
    monorepo_root is None,
    reason=(
        "Requires the full cardiacFoam monorepo tree (tutorials/). "
        "Clone the full repository to enable verification contracts tests."
    ),
)
