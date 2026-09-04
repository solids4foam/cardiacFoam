"""Behavioral test for ionicHeterogeneity::synthesizeTransmuralBandRegions.

This is a pure, mesh-free utility function in src/ionicModels/ionicModel/
ionicHeterogeneity.{H,C} used to give eikonalECG (a downstream consumer) the
same 3-named-region view of transmuralBands that
ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity already
builds inline for the monodomain path (ionicHeterogeneityOrchestrator.C,
"transmuralBands is a shorthand: expand to the same three named regions").
This test does not touch or depend on that orchestrator file; it only
verifies the new standalone function and cross-checks it against the
already-existing transmuralBandWeights() function, proving the two weighting
paths agree to numerical precision -- the "N=3 equivalence" property the
namedRegions generalization plan relies on.

Compiles a small throwaway C++ program against the already-built
libionicModels/libOpenFOAM shared libraries (same approach as
src/electroModels/tests/test_eikonal_template_generator.py).

Requires an OpenFOAM v2412 environment to be sourced (WM_PROJECT_DIR,
FOAM_LIBBIN, FOAM_USER_LIBBIN) with libionicModels already built; skips
otherwise.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
HEADER = ROOT / "src/ionicModels/ionicModel/ionicHeterogeneity.H"
SOURCE = ROOT / "src/ionicModels/ionicModel/ionicHeterogeneity.C"
IONIC_MODEL_DIR = ROOT / "src/ionicModels/ionicModel"


_CPP_SOURCE = '''\
#include "ionicHeterogeneity.H"
#include "scalarList.H"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Foam;
using namespace Foam::ionicHeterogeneity;

int main()
{
    const scalar endoMInterface = 0.3;
    const scalar mEpiInterface = 0.7;
    const scalar transitionWidth = 0.1;
    const word smoothing = "smoothstep";
    const word transitionMode = "blend";

    const List<NamedFieldRegion> regions =
        synthesizeTransmuralBandRegions(endoMInterface, mEpiInterface);

    std::cout << regions.size() << std::endl;
    forAll(regions, i)
    {
        std::cout
            << regions[i].name << " "
            << std::setprecision(12) << regions[i].rangeMin << " "
            << regions[i].rangeMax << " "
            << regions[i].baseline << std::endl;
    }

    // Cross-check against transmuralBandWeights() at a sweep of t values,
    // including exact boundaries and blend-zone midpoints.
    const scalarList tValues
    {
        0.0, 0.15, 0.3, 0.35, 0.5, 0.7, 0.75, 0.85, 1.0
    };

    scalar maxDiff = 0.0;
    forAll(tValues, i)
    {
        const scalar t = tValues[i];

        const TransmuralBandWeights legacy = transmuralBandWeights
        (
            t, endoMInterface, mEpiInterface, transitionWidth, smoothing,
            transitionMode
        );

        const List<NamedRegionWeight> generalized = namedRegionWeightsAt
        (
            t, regions, transitionWidth, smoothing, transitionMode
        );

        scalar wEndo = 0.0, wMid = 0.0, wEpi = 0.0;
        forAll(generalized, j)
        {
            if (generalized[j].name == "endocardialCells")
            {
                wEndo = generalized[j].weight;
            }
            else if (generalized[j].name == "mCells")
            {
                wMid = generalized[j].weight;
            }
            else if (generalized[j].name == "epicardialCells")
            {
                wEpi = generalized[j].weight;
            }
        }

        maxDiff = max(maxDiff, mag(wEndo - legacy.endo));
        maxDiff = max(maxDiff, mag(wMid - legacy.mCell));
        maxDiff = max(maxDiff, mag(wEpi - legacy.epi));

        std::cout
            << "t=" << t
            << " legacy=(" << legacy.endo << "," << legacy.mCell << ","
            << legacy.epi << ")"
            << " generalized=(" << wEndo << "," << wMid << "," << wEpi << ")"
            << std::endl;
    }

    std::cout << "maxDiff " << std::setprecision(15) << maxDiff << std::endl;

    return 0;
}
'''


def _openfoam_build_env():
    wm_project_dir = os.environ.get("WM_PROJECT_DIR")
    foam_libbin = os.environ.get("FOAM_LIBBIN")
    foam_user_libbin = os.environ.get("FOAM_USER_LIBBIN")

    if not (wm_project_dir and foam_libbin and foam_user_libbin):
        return None

    wm_project_dir = Path(wm_project_dir)
    foam_libbin = Path(foam_libbin)
    foam_user_libbin = Path(foam_user_libbin)

    if not (wm_project_dir / "src/OpenFOAM/lnInclude").is_dir():
        return None

    if not list(foam_user_libbin.glob("libionicModels.*")):
        return None

    return {
        "wm_project_dir": wm_project_dir,
        "foam_libbin": foam_libbin,
        "foam_user_libbin": foam_user_libbin,
    }


def _compiler():
    for candidate in ("clang++", "c++"):
        found = shutil.which(candidate)
        if found:
            return found
    return None


def _build(tmp_path: Path, name: str) -> Path:
    env = _openfoam_build_env()
    if env is None:
        pytest.skip(
            "OpenFOAM v2412 environment not sourced or libionicModels not "
            "built"
        )

    compiler = _compiler()
    if compiler is None:
        pytest.skip("C++ compiler unavailable")

    cpp_file = tmp_path / f"{name}.cpp"
    obj_file = tmp_path / f"{name}.o"
    exe_file = tmp_path / name

    cpp_file.write_text(_CPP_SOURCE)

    compile_cmd = [
        compiler,
        "-std=c++17", "-m64", "-pthread", "-ftrapping-math",
        "-DOPENFOAM=2412", "-DWM_DP", "-DWM_LABEL_SIZE=32",
        "-O3", "-DNoRepository", "-ftemplate-depth-100",
        "-DOPENFOAM_COM", "-DOPENFOAM_NOT_EXTEND",
        "-I", str(IONIC_MODEL_DIR),
        "-I", str(env["wm_project_dir"] / "src/OpenFOAM/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/OSspecific/POSIX/lnInclude"),
        "-fPIC",
        "-c", str(cpp_file), "-o", str(obj_file),
    ]
    subprocess.run(compile_cmd, check=True, capture_output=True, text=True)

    link_cmd = [
        compiler,
        "-std=c++17", "-m64", "-pthread", "-O3",
        str(obj_file),
        "-L", str(env["foam_libbin"]),
        "-L", str(env["foam_user_libbin"]),
        "-Wl,-rpath," + str(env["foam_libbin"]),
        "-Wl,-rpath," + str(env["foam_user_libbin"]),
        "-lionicModels", "-lOpenFOAM",
        "-o", str(exe_file),
    ]
    subprocess.run(link_cmd, check=True, capture_output=True, text=True)

    return exe_file


def test_header_declares_the_new_function():
    text = HEADER.read_text()
    assert "synthesizeTransmuralBandRegions" in text


def test_synthesizes_three_regions_matching_orchestrator_construction(tmp_path):
    exe = _build(tmp_path, "synthesize_regions_ok")
    result = subprocess.run([str(exe)], check=True, capture_output=True, text=True)

    lines = result.stdout.strip().splitlines()
    assert lines[0] == "3"
    assert lines[1] == "endocardialCells 0 0.3 endocardialCells"
    assert lines[2] == "mCells 0.3 0.7 mCells"
    assert lines[3] == "epicardialCells 0.7 1 epicardialCells"


def test_weights_match_transmuralBandWeights_to_numerical_precision(tmp_path):
    exe = _build(tmp_path, "synthesize_regions_equivalence")
    result = subprocess.run([str(exe)], check=True, capture_output=True, text=True)

    last_line = result.stdout.strip().splitlines()[-1]
    assert last_line.startswith("maxDiff")
    max_diff = float(last_line.split()[1])
    assert max_diff < 1e-12, result.stdout
