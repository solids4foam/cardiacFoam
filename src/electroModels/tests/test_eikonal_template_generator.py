"""Behavioral test for eikonalTemplateGenerator::generatePersonalizedTemplates.

This is a standalone component (src/electroModels/ecgModels/eikonalECG/
eikonalTemplateGenerator.{H,C}) that is not wired into eikonalECG yet (that is
a separate, later task). Because the function under test constructs a real
Foam::ionicModel (dictionary, autoPtr, scalarField, ...), it cannot be
exercised with a plain throwaway g++ program the way a header-only CellML
model can (see test_land_niederer_batched_preconditioning.py) -- it needs to
be compiled and linked against the already-built libelectroModels/
libionicModels/libOpenFOAM shared libraries.

Rather than adding a permanent applications/utilities/* wmake application
(this repo's convention for OpenFOAM-linked executables, see
applications/utilities/ionicHeterogeneityProbe), this test compiles a small
throwaway C++ program directly with the same compiler/linker flags `wmake`
used to build src/electroModels (captured from a real build), runs it, and
asserts on its stdout. This keeps the test self-contained in this one file.

Requires an OpenFOAM v2412 environment to be sourced (WM_PROJECT_DIR,
FOAM_LIBBIN, FOAM_USER_LIBBIN) with libelectroModels/libionicModels already
built (see AGENT_GUIDE.md / cardiacFoam build docs); skips otherwise.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
HEADER = ROOT / "src/electroModels/ecgModels/eikonalECG/eikonalTemplateGenerator.H"
SOURCE = ROOT / "src/electroModels/ecgModels/eikonalECG/eikonalTemplateGenerator.C"
MAKE_FILES = ROOT / "src/electroModels/Make/files"
ECG_INCLUDE_DIR = ROOT / "src/electroModels/ecgModels/eikonalECG"


# A minimal, physiologically-sane ionicModelConfig/heterogeneityDict pair.
# BuenoOroviocompactBatched is the runtime-selectable name for the batched
# Bueno-Orovio model (BuenoOrovioBatched itself is an unregistered base
# class -- only its "compactBatched" variant is added to the ionicModel
# runtime-selection table); it is one of the four *Batched families that
# support transmuralBands heterogeneity, and its 4-state ODE is cheap.
_IONIC_MODEL_CONFIG = """\
ionicModel        BuenoOroviocompactBatched;
tissue            epicardialCells;
batchedIntegrator rushLarsen;
batchedSubsteps   10;
singleCellStimulus
{
    stim_start      20;
    stim_duration   1;
    stim_amplitude  0.4;
    stim_period_S1  1000;
    nstim1          1;
    stim_period_S2  0;
    nstim2          0;
}
"""

_HETEROGENEITY_DICT = """\
mode             transmuralBands;
endoMInterface   0.3;
mEpiInterface    0.7;
transitionWidth  0.1;
transitionMode   blend;
smoothing        smoothstep;
"""

# Regression fixture for Finding 1 (mid-myocardium reference point landing in
# the blend zone). With endoMInterface=0.3, mEpiInterface=0.7, the pure
# mid-myocardium region under blend mode is [endoMInterface+transitionWidth,
# mEpiInterface]. transitionWidth=0.25 -> pure region [0.55, 0.7]; the old
# buggy tPoints[1]=0.5*(endoM+mEpi)=0.5 formula landed inside the endo/mid
# blend zone [0.3, 0.55) instead. transitionWidth=0.05 -> pure region
# [0.35, 0.7], which the buggy midpoint 0.5 already happened to fall inside,
# so under the OLD code both widths gave the same weighting for the "mid"
# anchor; under the FIXED code both widths now correctly resolve to a point
# in the pure mid-myocardium region, so their "mid" templates should match.
_HETEROGENEITY_DICT_WIDE_TRANSITION = """\
mode             transmuralBands;
endoMInterface   0.3;
mEpiInterface    0.7;
transitionWidth  0.25;
transitionMode   blend;
smoothing        smoothstep;
"""

_HETEROGENEITY_DICT_NARROW_TRANSITION = """\
mode             transmuralBands;
endoMInterface   0.3;
mEpiInterface    0.7;
transitionWidth  0.05;
transitionMode   blend;
smoothing        smoothstep;
"""

# Fatal-path fixture for Finding 1's new guard: endoMInterface=0.3,
# mEpiInterface=0.5, transitionWidth=0.2 makes endoMInterface+transitionWidth
# exactly equal to mEpiInterface (0.3+0.2==0.5 is exact in IEEE754 double),
# i.e. a zero-width pure mid-myocardium region. This is deliberately NOT the
# same config that ionicHeterogeneity::validateTransmuralBandConfig already
# rejects (that check is endoMInterface+transitionWidth > mEpiInterface+SMALL,
# i.e. it only fires when the transition bands *overlap*; it does not fire on
# a zero-width -- but non-overlapping -- pure region, which is the gap this
# new guard closes). Verified empirically: endoMInterface=0.3,
# mEpiInterface=0.7, transitionWidth=0.45 (as a naive "clearly too wide"
# choice) actually still trips the pre-existing overlap check first, not the
# new guard -- so this fixture was chosen specifically to isolate the new
# guard's fatal path.
_HETEROGENEITY_DICT_ZERO_WIDTH_MID = """\
mode             transmuralBands;
endoMInterface   0.3;
mEpiInterface    0.5;
transitionWidth  0.2;
transitionMode   blend;
smoothing        smoothstep;
"""


def _cpp_source(capture_duration: str, heterogeneity_dict: str = _HETEROGENEITY_DICT) -> str:
    return f'''\
#include "eikonalTemplateGenerator.H"
#include "IStringStream.H"
#include "dictionary.H"

#include <iomanip>
#include <iostream>

using namespace Foam;
using namespace Foam::eikonalECG_templates;

namespace
{{

dictionary parseDict(const std::string& text)
{{
    IStringStream is(text);
    return dictionary(is);
}}

void report(const word& name, const DynamicTemplate& tpl)
{{
    const label n = tpl.times.size();
    scalar resting = tpl.valuesMv[0];
    scalar peak = resting;
    forAll(tpl.valuesMv, i)
    {{
        peak = max(peak, tpl.valuesMv[i]);
    }}
    std::cout
        << name << " " << n << " "
        << std::setprecision(12) << tpl.times[0] << " "
        << tpl.times[n - 1] << " "
        << resting << " " << peak << " "
        << tpl.valuesMv[n - 1] << std::endl;
}}

}} // End unnamed namespace


int main()
{{
    const std::string ionicModelConfigText =
{_quote_cpp_lines(_IONIC_MODEL_CONFIG)};

    const std::string heterogeneityDictText =
{_quote_cpp_lines(heterogeneity_dict)};

    const dictionary ionicModelConfig(parseDict(ionicModelConfigText));
    const dictionary heterogeneityDict(parseDict(heterogeneityDictText));

    const label nBeats = 2;
    const scalar captureDuration = {capture_duration};
    const scalar dt = 0.0001;

    const TemplateTriplet templates =
        generatePersonalizedTemplates
        (
            ionicModelConfig,
            heterogeneityDict,
            nBeats,
            captureDuration,
            dt
        );

    report("endo", templates.endo);
    report("mid", templates.mid);
    report("epi", templates.epi);

    return 0;
}}
'''


def _quote_cpp_lines(text: str) -> str:
    lines = text.splitlines(keepends=True)
    escaped = [line.replace("\\", "\\\\").replace('"', '\\"') for line in lines]
    escaped = [line.replace("\n", "\\n") for line in escaped]
    return "\n".join(f'    "{line}"' for line in escaped)


def _openfoam_build_env():
    """Return the OpenFOAM env vars/paths needed to compile+link, or None."""
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

    required_libs = ["libelectroModels", "libionicModels"]
    for lib in required_libs:
        candidates = list(foam_user_libbin.glob(f"{lib}.*"))
        if not candidates:
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


def _build(
    tmp_path: Path,
    name: str,
    capture_duration: str,
    heterogeneity_dict: str = _HETEROGENEITY_DICT,
) -> Path:
    env = _openfoam_build_env()
    if env is None:
        pytest.skip(
            "OpenFOAM v2412 environment not sourced or "
            "libelectroModels/libionicModels not built"
        )

    compiler = _compiler()
    if compiler is None:
        pytest.skip("C++ compiler unavailable")

    cpp_file = tmp_path / f"{name}.cpp"
    obj_file = tmp_path / f"{name}.o"
    exe_file = tmp_path / name

    cpp_file.write_text(_cpp_source(capture_duration, heterogeneity_dict))

    # These flags mirror src/electroModels/Make/options (and
    # src/ionicModels/Make/options for library-specific flags such as
    # -DWM_LABEL_SIZE=32 and the -l link list below). If those change,
    # re-sync this list by hand -- there is no automated derivation here.
    compile_cmd = [
        compiler,
        "-std=c++17", "-m64", "-pthread", "-ftrapping-math",
        "-DOPENFOAM=2412", "-DWM_DP", "-DWM_LABEL_SIZE=32",
        "-O3", "-DNoRepository", "-ftemplate-depth-100",
        "-DOPENFOAM_COM", "-DOPENFOAM_NOT_EXTEND",
        "-I", str(ECG_INCLUDE_DIR),
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
        "-lelectroModels", "-lionicModels", "-lgenericWriter",
        "-lactiveTensionModels", "-lphysicsModel",
        "-lfiniteVolume", "-lmeshTools", "-lsurfMesh",
        "-ldynamicFvMesh", "-ldynamicMesh", "-lODE", "-lOpenFOAM",
        "-o", str(exe_file),
    ]
    subprocess.run(link_cmd, check=True, capture_output=True, text=True)

    return exe_file


def test_header_exposes_only_the_documented_api():
    text = HEADER.read_text()
    assert "namespace eikonalECG_templates" in text
    assert "struct DynamicTemplate" in text
    assert "struct TemplateTriplet" in text
    assert "generatePersonalizedTemplates" in text


def test_generator_does_not_touch_eikonalECG_or_tissueTemplates():
    source = SOURCE.read_text()
    header = HEADER.read_text()
    for forbidden in ('"eikonalECG.H"', '"tissueTemplates.H"'):
        assert forbidden not in source
        assert forbidden not in header


def test_make_files_lists_the_new_source():
    text = MAKE_FILES.read_text()
    assert "ecgModels/eikonalECG/eikonalECG.C" in text
    assert "ecgModels/eikonalECG/eikonalTemplateGenerator.C" in text


def test_generates_three_distinct_valid_templates(tmp_path):
    exe = _build(tmp_path, "eikonal_template_generator_ok", "0.6")

    result = subprocess.run(
        [str(exe)], check=True, capture_output=True, text=True
    )

    lines = [line.split() for line in result.stdout.strip().splitlines()]
    assert len(lines) == 3

    traces = {}
    for name, n, t0, tLast, resting, peak, final in lines:
        n = int(n)
        t0, tLast, resting, peak, final = map(
            float, (t0, tLast, resting, peak, final)
        )
        traces[name] = dict(
            n=n, t0=t0, tLast=tLast, resting=resting, peak=peak, final=final
        )

    assert set(traces) == {"endo", "mid", "epi"}

    for name, tr in traces.items():
        # ceil(0.6/0.0001) + 1
        assert tr["n"] == 6001, name
        assert tr["t0"] == pytest.approx(0.0, abs=1e-9), name
        assert tr["tLast"] == pytest.approx(0.6, abs=1e-9), name
        # Physiological resting range and a comfortable upstroke amplitude.
        assert -95.0 < tr["resting"] < -75.0, name
        assert tr["peak"] - tr["resting"] >= 50.0, name
        # Recovered to within the 5 mV validation tolerance (with margin).
        assert abs(tr["final"] - tr["resting"]) < 5.0, name

    # Heterogeneity actually took effect: the three anchors are not
    # identical copies of one another.
    peaks = [traces[k]["peak"] for k in ("endo", "mid", "epi")]
    assert len(set(round(p, 3) for p in peaks)) == 3


def test_fatals_naming_the_anchor_when_capture_is_too_short(tmp_path):
    exe = _build(tmp_path, "eikonal_template_generator_short", "0.05")

    result = subprocess.run([str(exe)], capture_output=True, text=True)

    assert result.returncode != 0
    assert "did not recover to within 5 mV" in result.stderr
    assert any(
        anchor in result.stderr
        for anchor in ("endocardium", "mid-myocardium", "epicardium")
    )


def _mid_peak(exe: Path) -> float:
    result = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
    for line in result.stdout.strip().splitlines():
        fields = line.split()
        if fields[0] == "mid":
            # fields: name n t0 tLast resting peak final
            return float(fields[5])
    raise AssertionError(f"no 'mid' line in stdout: {result.stdout!r}")


def test_wide_transition_width_still_yields_a_pure_mid_myocardium_template(
    tmp_path,
):
    # Regression test for Finding 1: eikonalTemplateGenerator.C used to place
    # the mid-myocardium reference point at the plain midpoint
    # 0.5*(endoMInterface+mEpiInterface), which can land inside the
    # endo<->mid blend zone instead of the actual pure mid-myocardium region
    # once transitionWidth is wide enough (see the
    # _HETEROGENEITY_DICT_WIDE_TRANSITION / _NARROW_TRANSITION docstrings
    # above for the exact interval arithmetic).
    #
    # Under the OLD buggy code, tPoints[1] was always exactly 0.5 regardless
    # of transitionWidth, so widening transitionWidth from 0.05 to 0.25 moved
    # that fixed point from the pure mid region into the blend zone, and the
    # resulting "mid" anchor's simulated peak Vm would shift accordingly.
    # Under the FIXED code, both widths resolve tPoints[1] into the pure mid
    # region ([0.35, 0.7] and [0.55, 0.7] respectively), which the ionic
    # model resolves to the identical mid-myocardium parameterization
    # (transmuralBandWeights returns {0, 1, 0} anywhere in that band), so the
    # two "mid" templates' peak Vm should match to numerical precision.
    exe_wide = _build(
        tmp_path,
        "eikonal_template_generator_wide_transition",
        "0.6",
        _HETEROGENEITY_DICT_WIDE_TRANSITION,
    )
    exe_narrow = _build(
        tmp_path,
        "eikonal_template_generator_narrow_transition",
        "0.6",
        _HETEROGENEITY_DICT_NARROW_TRANSITION,
    )

    peak_wide = _mid_peak(exe_wide)
    peak_narrow = _mid_peak(exe_narrow)

    assert peak_wide == pytest.approx(peak_narrow, abs=1e-6)


def test_fatals_when_transition_width_leaves_no_pure_mid_region(tmp_path):
    # Regression test for Finding 1's new validation guard: endoMInterface=
    # 0.3, mEpiInterface=0.5, transitionWidth=0.2 makes
    # endoMInterface+transitionWidth exactly equal to mEpiInterface, i.e. a
    # zero-width pure mid-myocardium region -- a configuration that
    # ionicHeterogeneity::validateTransmuralBandConfig's overlap check does
    # NOT catch (that check only fires once the transition bands actually
    # overlap), so the new guard in eikonalTemplateGenerator.C is what must
    # catch it.
    exe = _build(
        tmp_path,
        "eikonal_template_generator_zero_width_mid",
        "0.6",
        _HETEROGENEITY_DICT_ZERO_WIDTH_MID,
    )

    result = subprocess.run([str(exe)], capture_output=True, text=True)

    assert result.returncode != 0
    assert "transitionWidth" in result.stderr
    assert "pure mid-myocardium region" in result.stderr
