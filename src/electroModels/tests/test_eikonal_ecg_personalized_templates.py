"""Behavioral test for eikonalECG's personalizedTemplates wiring (Task 2).

Task 1 (test_eikonal_template_generator.py) built and verified the standalone
Foam::eikonalECG_templates::generatePersonalizedTemplates(...) component in
isolation. This task wires that component into Foam::eikonalECG: a new
"personalizedTemplates" sub-dictionary on the eikonalECG domain dict, the
constructor-level "reject before construction" validation the plan's shared
Configuration contract requires, generatePersonalizedTemplates() (the member
function, not to be confused with the Task 1 free function it calls), and a
new branch in reconstructGradVm() that evaluates the generated templates
instead of the compiled-in tissueTemplates.H arrays.

eikonalECG's constructor takes only a dictionary& (no mesh/Time/ecgDomain), so
every "reject before construction" check that does not depend on the case's
runtime ionicHeterogeneity field (i.e. everything except "missing
ionicHeterogeneity block" / "mode != transmuralBands", which are only
resolvable once a real mesh + electroProperties IOdictionary exist) can be
exercised end-to-end here by constructing Foam::eikonalECG directly from a
hand-built dictionary and observing whether it fatals. This test compiles a
small throwaway C++ program against the already-built libelectroModels (same
approach and captured wmake flags as test_eikonal_template_generator.py) that
does exactly that.

What this file does NOT cover (needs a real fvMesh/ecgDomain, deferred to
Task 3's tutorial-based regression coverage per the task brief):
  - generatePersonalizedTemplates() member's own mesh-dependent checks
    (missing ionicHeterogeneity block; mode != transmuralBands), which read
    constant/electroProperties off a real mesh.
  - The new dynamic-template branch in reconstructGradVm() actually producing
    correct field values end-to-end (verified only by source inspection here
    that the branch exists, selects the generated triplet, and applies the
    same mV->V unit conversion as the compiled-array branch).
  - solve()'s call-once behavior (personalizedTemplatesGenerated_ guard).

Requires an OpenFOAM v2412 environment to be sourced (WM_PROJECT_DIR,
FOAM_LIBBIN, FOAM_USER_LIBBIN) with libelectroModels already built; skips
otherwise.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
HEADER = ROOT / "src/electroModels/ecgModels/eikonalECG/eikonalECG.H"
SOURCE = ROOT / "src/electroModels/ecgModels/eikonalECG/eikonalECG.C"
ELECTRO_MODELS = ROOT / "src/electroModels"

TASK1_FILES = [
    ELECTRO_MODELS / "ecgModels/eikonalECG/eikonalTemplateGenerator.H",
    ELECTRO_MODELS / "ecgModels/eikonalECG/eikonalTemplateGenerator.C",
    ELECTRO_MODELS / "ecgModels/eikonalECG/tissueTemplates.H",
]


# A minimal, valid singleCellStimulus + ionicModelConfig, matching the
# Configuration contract's example (Task 1's BuenoOroviocompactBatched
# choice: the runtime-selectable "compactBatched" variant of the batched
# Bueno-Orovio model; cheap 4-state ODE, one of the four *Batched families
# supporting transmuralBands heterogeneity).
_VALID_IONIC_MODEL_CONFIG = """\
ionicModelConfig
{
    ionicModel BuenoOroviocompactBatched;
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
}
"""


def _domain_dict(
    personalized_templates_body: str = None,
    extra_top_level: str = "",
) -> str:
    """Build a minimal eikonalECG domain dict (the constructor's argument).

    personalized_templates_body, if given, is embedded verbatim as the
    "personalizedTemplates { ... }" sub-dict's contents; if None, the
    sub-dict is omitted entirely (personalizedTemplates disabled).
    """
    block = ""
    if personalized_templates_body is not None:
        block = f"personalizedTemplates\n{{\n{personalized_templates_body}\n}}\n"

    return f"""\
sampling
{{
    start   0;
    end     0.1;
    deltaT  0.001;
}}
{extra_top_level}
{block}
"""


_VALID_PERSONALIZED_TEMPLATES = f"""\
{_VALID_IONIC_MODEL_CONFIG}
nBeats   2;
duration 0.6;
dt       0.0001;
"""


def _cpp_source(dict_text: str) -> str:
    return f'''\
#include "eikonalECG.H"
#include "IStringStream.H"
#include "dictionary.H"

#include <iostream>

using namespace Foam;

namespace
{{

dictionary parseDict(const std::string& text)
{{
    IStringStream is(text);
    return dictionary(is);
}}

}} // End unnamed namespace


int main()
{{
    const std::string dictText =
{_quote_cpp_lines(dict_text)};

    const dictionary dict(parseDict(dictText));
    eikonalECG solver(dict);

    std::cout << "constructed_ok" << std::endl;
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

    if not (ELECTRO_MODELS / "lnInclude" / "eikonalECG.H").exists():
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


def _build(tmp_path: Path, name: str, dict_text: str) -> Path:
    env = _openfoam_build_env()
    if env is None:
        pytest.skip(
            "OpenFOAM v2412 environment not sourced or libelectroModels not "
            "built (with an up-to-date lnInclude)"
        )

    compiler = _compiler()
    if compiler is None:
        pytest.skip("C++ compiler unavailable")

    cpp_file = tmp_path / f"{name}.cpp"
    obj_file = tmp_path / f"{name}.o"
    exe_file = tmp_path / name

    cpp_file.write_text(_cpp_source(dict_text))

    # These flags mirror src/electroModels/Make/options; the electroModels
    # package's own generated lnInclude/ (a flat directory of symlinks to
    # every header the library compiles, including eikonalECG.H) stands in
    # for the many individual -I<subdir> entries wmake would otherwise pass,
    # so only the cross-package OpenFOAM lnInclude dirs need listing by hand.
    # If Make/options changes, re-sync this list.
    compile_cmd = [
        compiler,
        "-std=c++17", "-m64", "-pthread", "-ftrapping-math",
        "-DOPENFOAM=2412", "-DWM_DP", "-DWM_LABEL_SIZE=32",
        "-O3", "-DNoRepository", "-ftemplate-depth-100",
        "-DOPENFOAM_COM", "-DOPENFOAM_NOT_EXTEND",
        "-Wno-undefined-var-template",
        "-I", str(ELECTRO_MODELS / "lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/finiteVolume/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/meshTools/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/surfMesh/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/dynamicMesh/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/dynamicFvMesh/lnInclude"),
        "-I", str(env["wm_project_dir"] / "src/ODE/lnInclude"),
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


def _run(tmp_path, name, dict_text):
    exe = _build(tmp_path, name, dict_text)
    return subprocess.run([str(exe)], capture_output=True, text=True)


# --------------------------------------------------------------------------
# Static checks
# --------------------------------------------------------------------------


def test_header_declares_personalized_template_members():
    text = HEADER.read_text()
    assert "personalizedTemplatesEnabled_" in text
    assert "personalizedTemplatesGenerated_" in text
    assert "personalizedTemplatesDict_" in text
    assert "eikonalECG_templates::TemplateTriplet personalizedTemplates_" in text
    assert "void generatePersonalizedTemplates(const ecgDomain& domain);" in text
    assert '#include "eikonalTemplateGenerator.H"' in text


def test_header_declares_shared_heterogeneity_lookup_helper():
    text = HEADER.read_text()
    assert "findHeterogeneityDict" in text


def test_source_does_not_duplicate_heterogeneity_dict_lookup():
    # calculateTransmuralWeights() must delegate to findHeterogeneityDict()
    # rather than keep its own inline copy of the lookup, so it can never
    # silently disagree with generatePersonalizedTemplates() about which
    # ionicHeterogeneity dictionary is authoritative for a given case.
    text = SOURCE.read_text()
    assert text.count('electroProperties.findDict("ionicHeterogeneity")') == 1
    assert "findHeterogeneityDict(electroProperties)" in text


def test_reconstruct_grad_vm_has_a_personalized_template_branch():
    text = SOURCE.read_text()
    assert "personalizedTemplatesEnabled_" in text
    assert "personalizedTemplates_.endo" in text
    assert "personalizedTemplates_.mid" in text
    assert "personalizedTemplates_.epi" in text
    # Same mV->V unit convention as the pre-existing compiled-array branch.
    assert text.count("rawDUds*1e-3") + text.count("rawDUds * 1e-3") == 2


def test_solve_generates_templates_once_when_opted_in():
    text = SOURCE.read_text()
    assert (
        "personalizedTemplatesEnabled_ && !personalizedTemplatesGenerated_"
        in text
    )


def test_task1_files_untouched():
    # This task must not modify eikonalTemplateGenerator.{H,C} or
    # tissueTemplates.H; those were reviewed and closed under Task 1.
    result = subprocess.run(
        ["git", "status", "--porcelain", "--"] + [str(p) for p in TASK1_FILES],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == "", result.stdout


# --------------------------------------------------------------------------
# Dynamic checks: eikonalECG(dict) construction (no mesh required)
# --------------------------------------------------------------------------


def test_valid_personalized_templates_config_constructs(tmp_path):
    dict_text = _domain_dict(_VALID_PERSONALIZED_TEMPLATES)
    result = _run(tmp_path, "valid_config", dict_text)
    assert result.returncode == 0, result.stderr
    assert "constructed_ok" in result.stdout


def test_no_personalized_templates_block_still_constructs(tmp_path):
    dict_text = _domain_dict(None)
    result = _run(tmp_path, "no_block", dict_text)
    assert result.returncode == 0, result.stderr
    assert "constructed_ok" in result.stdout


def test_rejects_missing_ionic_model_config(tmp_path):
    body = "nBeats   2;\nduration 0.6;\ndt       0.0001;\n"
    result = _run(tmp_path, "missing_ionic_model_config", _domain_dict(body))
    assert result.returncode != 0
    assert "ionicModelConfig" in result.stderr


def test_rejects_ionic_model_config_missing_ionic_model_key(tmp_path):
    body = f"""\
ionicModelConfig
{{
    singleCellStimulus
    {{
        stim_start      20;
        stim_duration   1;
        stim_amplitude  0.4;
        stim_period_S1  1000;
        nstim1          1;
        stim_period_S2  0;
        nstim2          0;
    }}
}}
nBeats   2;
duration 0.6;
dt       0.0001;
"""
    result = _run(tmp_path, "missing_ionic_model_key", _domain_dict(body))
    assert result.returncode != 0
    assert "'ionicModel' entry" in result.stderr


def test_rejects_ionic_model_config_missing_single_cell_stimulus(tmp_path):
    body = """\
ionicModelConfig
{
    ionicModel BuenoOroviocompactBatched;
}
nBeats   2;
duration 0.6;
dt       0.0001;
"""
    result = _run(tmp_path, "missing_single_cell_stimulus", _domain_dict(body))
    assert result.returncode != 0
    assert "singleCellStimulus" in result.stderr


def test_rejects_nbeats_below_one(tmp_path):
    body = f"""\
{_VALID_IONIC_MODEL_CONFIG}
nBeats   0;
duration 0.6;
dt       0.0001;
"""
    result = _run(tmp_path, "nbeats_zero", _domain_dict(body))
    assert result.returncode != 0
    assert "nBeats" in result.stderr


@pytest.mark.parametrize("duration,dt", [("0", "0.0001"), ("0.6", "0")])
def test_rejects_non_positive_duration_or_dt(tmp_path, duration, dt):
    body = f"""\
{_VALID_IONIC_MODEL_CONFIG}
nBeats   2;
duration {duration};
dt       {dt};
"""
    result = _run(tmp_path, f"non_positive_{duration}_{dt}", _domain_dict(body))
    assert result.returncode != 0
    assert "positive" in result.stderr


def test_rejects_duration_exceeding_one_s1_period(tmp_path):
    # stim_period_S1=1000 ms => one period is 1.0 s; duration=2.0 s exceeds
    # it, so capture would no longer contain exactly the final S1 response.
    body = f"""\
{_VALID_IONIC_MODEL_CONFIG}
nBeats   2;
duration 2.0;
dt       0.0001;
"""
    result = _run(tmp_path, "duration_exceeds_period", _domain_dict(body))
    assert result.returncode != 0
    assert "S1 period" in result.stderr


def test_rejects_s2_pacing(tmp_path):
    body = """\
ionicModelConfig
{
    ionicModel BuenoOroviocompactBatched;
    singleCellStimulus
    {
        stim_start      20;
        stim_duration   1;
        stim_amplitude  0.4;
        stim_period_S1  1000;
        nstim1          1;
        stim_period_S2  500;
        nstim2          1;
    }
}
nBeats   2;
duration 0.6;
dt       0.0001;
"""
    result = _run(tmp_path, "s2_pacing", _domain_dict(body))
    assert result.returncode != 0
    assert "S2 pacing" in result.stderr


def test_rejects_manufactured_ecg_combined_with_personalized_templates(
    tmp_path,
):
    dict_text = _domain_dict(
        _VALID_PERSONALIZED_TEMPLATES,
        extra_top_level="manufacturedEikonalECG {}\n",
    )
    result = _run(tmp_path, "manufactured_conflict", dict_text)
    assert result.returncode != 0
    assert "manufactured ECG" in result.stderr


def test_rejects_manufactured_ecg_verifier_combined_with_personalized_templates(
    tmp_path,
):
    dict_text = _domain_dict(
        _VALID_PERSONALIZED_TEMPLATES,
        extra_top_level=(
            "verificationModel\n"
            "{\n"
            "    type manufacturedEikonalECGVerifier;\n"
            "}\n"
        ),
    )
    result = _run(tmp_path, "manufactured_verifier_conflict", dict_text)
    assert result.returncode != 0
    assert "manufactured ECG" in result.stderr
