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
#     test_runtime_dependencies_live_verification
#
# Description
#     Cross-checks resolve_runtime_dependencies against a real sourced
#     OpenFOAM install, confirming the I3b findings the resolver is built on
#     rather than assuming them. Skipped when OpenFOAM is not sourced.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""The live ratchet for I3b: does the resolver actually behave the way the
plan's verified-live findings say a sourced OpenFOAM install behaves?

To run it::

    source /Volumes/OpenFOAM-v2412/etc/bashrc
    cd applications/scripts/driverFoam
    uv run pytest openfoam_driver/tests/plugins/cardiacfoam/\\
test_runtime_dependencies_live_verification.py -v

Skips (not fails) when ``cardiacFoam`` is not on ``PATH`` -- the normal state
without a sourced OpenFOAM install. A skip is not a pass.
"""

from __future__ import annotations

import os
import shutil

import pytest

from openfoam_driver.plugins.cardiacfoam.runtime_evidence import (
    resolve_runtime_dependencies,
)

requires_openfoam = pytest.mark.skipif(
    shutil.which("cardiacFoam") is None,
    reason=(
        "cardiacFoam not on PATH -- source the OpenFOAM bashrc to run the "
        "live runtime-dependency check. This is a SKIP, not a pass."
    ),
)


@requires_openfoam
def test_foam_module_libbin_is_unset_in_this_sourced_shell() -> None:
    """I3b finding 1. If this ever starts failing, the resolver's core
    justification -- never map a name to one fixed lib variable -- needs
    re-examining, not just this test."""
    assert "FOAM_MODULE_LIBBIN" not in os.environ
    assert "FOAM_USER_LIBBIN" in os.environ


@requires_openfoam
def test_the_solver_and_its_always_linked_libraries_resolve(tmp_path) -> None:
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case")}

    assert deps["cardiacFoam"].path is not None
    assert deps["cardiacFoam"].path.is_file()

    for name in ("electroModels", "ionicModels", "genericWriter", "activeTensionModels"):
        assert deps[name].path is not None, f"{name} did not resolve"
        assert deps[name].path.suffix == ".dylib"


@requires_openfoam
def test_physicsmodel_resolves_despite_make_files_declaring_the_other_variable(
    tmp_path,
) -> None:
    """I3b finding 2: Make/files declares FOAM_MODULE_LIBBIN, but the
    installed library actually lives in FOAM_USER_LIBBIN."""
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case")}
    assert deps["physicsModel"].path is not None
    user_libbin = os.environ.get("FOAM_USER_LIBBIN", "")
    assert str(deps["physicsModel"].path).startswith(user_libbin)


@requires_openfoam
def test_electromechanicalmodels_is_genuinely_absent_in_this_install(tmp_path) -> None:
    """I3b finding 4: absent entirely in the maintainer's lightweight
    default -- reported unavailable, not an error, and not omitted."""
    deps = {d.name: d for d in resolve_runtime_dependencies(tmp_path / "case")}
    assert "electroMechanicalModels" in deps
    assert deps["electroMechanicalModels"].path is None
    assert deps["electroMechanicalModels"].required is False


@requires_openfoam
def test_a_real_manufactured_solution_tutorials_controldict_resolves(tmp_path) -> None:
    """I2c against the real repository file, not a synthetic fixture:
    monodomainTotalLagrangianEM's controlDict declares libverificationModels
    (bare) and a case-local .so built from sources inside the case. The
    case-local .so is not built in this checkout, so it must resolve to
    unavailable rather than a guessed path -- never a crash."""
    from openfoam_driver.tests.conftest import monorepo_root

    if monorepo_root is None:
        pytest.skip("requires the full cardiacFoam monorepo tree")

    case_root = (
        monorepo_root
        / "tutorials"
        / "manufacturedSolutions"
        / "monodomainTotalLagrangianEM"
    )
    if not (case_root / "system" / "controlDict").is_file():
        pytest.skip("monodomainTotalLagrangianEM tutorial not present in this checkout")

    deps = {d.name: d for d in resolve_runtime_dependencies(case_root)}

    assert deps["verificationModels"].required is True
    assert deps["verificationModels"].path is not None
    assert deps["verificationModels"].path.is_file()

    assert "manufacturedMonodomainTotalLagrangianEM" in deps
    case_local = deps["manufacturedMonodomainTotalLagrangianEM"]
    assert case_local.required is True
    # Not built in this checkout (no platforms/ dir under the tutorial) --
    # must be unavailable, never a fabricated path.
    if not (case_root / "platforms").exists():
        assert case_local.path is None
    else:
        assert case_local.path is None or case_local.path.is_file()
