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
#     test_mesh_geometry_contract
#
# Description
#     Guard against drift between the Python detector thresholds and the
#     checkMeshGeometry.C source-of-truth constants.
#
#     Parses the *named* C++ constants (mmLower/umLower/umUpper) rather than
#     matching bare literals, which could match unrelated numbers anywhere in
#     the source.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from pathlib import Path

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.core.specs import mesh_geometry


def _checkmesh_source() -> str:
    # tests/ -> openfoam_driver/ -> driverFoam/ -> scripts/ -> applications/
    applications = Path(__file__).resolve().parents[5]
    src = applications / "utilities" / "checkMeshGeometry" / "checkMeshGeometry.C"
    return src.read_text()


def _cpp_scalar(source: str, name: str) -> float:
    """Extract `const scalar <name> = <value>;` from the C++ source."""
    m = re.search(rf"\bscalar\s+{name}\s*=\s*([0-9.eE+\-]+)\s*;", source)
    if not m:
        pytest.fail(f"named C++ constant '{name}' not found in checkMeshGeometry.C")
    return float(m.group(1))


def test_named_cpp_constants_match_python():
    source = _checkmesh_source()
    assert _cpp_scalar(source, "mmLower") == mesh_geometry._MM_LOWER == 20.0
    assert _cpp_scalar(source, "umLower") == mesh_geometry._UM_LOWER == 1000.0
    assert _cpp_scalar(source, "umUpper") == mesh_geometry._UM_UPPER == 1.0e6


def test_python_scale_factors():
    assert mesh_geometry.classify_scale(50.0).scale_factor == 1e-3      # mm
    assert mesh_geometry.classify_scale(5000.0).scale_factor == 1e-6    # um
    assert mesh_geometry.classify_scale(0.1).scale_factor == 1.0        # m
    assert mesh_geometry.classify_scale(1.7).scale_factor == 1.0        # torso m
