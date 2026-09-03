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
#     test_eikonal_ecg_personalized_templates_regression
#
# Description
#     End-to-end driverFOAM regression coverage for eikonalECG's opt-in
#     personalizedTemplates generator (Task 3b of the personalized
#     action-potential template plan; see
#     .superpowers/sdd/task-3-brief.md). Five checks, driven entirely
#     through `driverFoam plan --strict` + `driverFoam run --run-document`
#     against copied temporary cases:
#
#       1. Fallback parity   -- the pre-existing manufacturedEikonalECG
#                                regression (no personalizedTemplates block)
#                                still reproduces its committed reference.
#       2. Opt-in smoke test -- the eikonalECGPersonalized tutorial (as
#                                committed) completes, writes a finite
#                                postProcessing/eikonalECG.dat, and generates
#                                exactly three action-potential templates.
#       3. Baseline parity   -- personalizedTemplates-present vs. the
#                                compiled-fallback (block absent) agree
#                                closely away from the upstroke transient.
#       4. Parameter sensitivity -- scaling a real TWorld constant
#                                (AC_GKr) through ionicConstantOverrides
#                                measurably changes the generated ECG.
#       5. Unit guard        -- the shared piecewise-linear derivative
#                                primitive, the mV/s->V/s conversion, and
#                                the chain-rule minus sign are each
#                                independently correct, tested
#                                deterministically (no ionic-model
#                                amplitude dependence).
#
#     Skips (not fails) when OpenFOAM/cardiacFoam is not sourced/built, or
#     when the full cardiacFoam monorepo tree is not present. A skip is not
#     a pass -- see the skip reasons for exactly what to source/build.
#
#     FIXED BUG (was "KNOWN, PRE-EXISTING BUG" as of this file's original
#     authorship, commits 32cd2436..08a514d9 on ep-work-onto-main; fixed on
#     top of 18d440d3): `driverFoam plan --strict` / `driverFoam run
#     --strict` used to incorrectly report ecgDomains.<name>.
#     personalizedTemplates.ionicModelConfig.ionicModel as "required" for
#     ANY ecgSolver=eikonalECG domain, even one with no personalizedTemplates
#     block at all -- because nBeats/duration/dt (but not ionicModel) carried
#     a catalog typical_value, so the dict-synthesis path's typical-value
#     fallback silently populated those three, making the fourth leaf's
#     required_when siblings look "present" even when the user never
#     configured personalizedTemplates. Fixed by removing typical_value from
#     nBeats/duration/dt in dict_entries_catalog.py (matching ionicModel,
#     which already had none) -- there is no case-independent default for
#     those three that should ever be silently synthesized into a case that
#     never opted into personalizedTemplates. `_plan_strict_tolerating_
#     known_fallback_bug` below (used by `fallback_result`) now requires a
#     clean "ok" plan like every other check in this file; the test that
#     used to pin this bug
#     (test_known_bug_plain_eikonal_ecg_domains_fail_strict_validation_due_
#     to_personalized_templates_catalog_defect) has been removed now that
#     the bug it pinned no longer reproduces.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


# --------------------------------------------------------------------------- #
# Environment gating (mirrors test_ionic_catalog_live_verification.py /
# test_runtime_dependencies_live_verification.py: skip, don't fail, when the
# OpenFOAM environment isn't sourced. A skip is not a pass.)
# --------------------------------------------------------------------------- #

requires_openfoam = pytest.mark.skipif(
    shutil.which("cardiacFoam") is None or not os.environ.get("WM_PROJECT_DIR"),
    reason=(
        "cardiacFoam not on PATH / WM_PROJECT_DIR unset -- source the "
        "OpenFOAM bashrc (and build libelectroModels) to run this live "
        "regression. This is a SKIP, not a pass."
    ),
)

pytestmark = [skip_without_monorepo, requires_openfoam]


def _compiler() -> str | None:
    for candidate in ("clang++", "c++"):
        found = shutil.which(candidate)
        if found:
            return found
    return None


# --------------------------------------------------------------------------- #
# driverFOAM subprocess helpers
# --------------------------------------------------------------------------- #

def _driver_env() -> dict[str, str]:
    """Environment for `python -m openfoam_driver ...` subprocess calls.

    Mirrors regression_equivalence/dual_run.py's `_drive_agent`, plus two
    backend-selection env vars this machine's build needs at runtime
    (documented in .superpowers/sdd/task-3a-report.md /
    task-3a-fix-report.md): FORCE_LIGHTWEIGHT_PHYSICSMODEL and
    DRIVERFOAM_CARDIACFOAM_BACKEND select the lightweight physicsModel
    backend that this checkout's cardiacFoam binary is actually linked
    against (verified via `otool -L` in task-3a-report.md), which
    `plan --strict`'s environment_preflight otherwise fails to detect.
    """
    import openfoam_driver

    pkg_parent = str(Path(openfoam_driver.__file__).resolve().parent.parent)
    env = os.environ.copy()
    env["PYTHONPATH"] = pkg_parent + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("FORCE_LIGHTWEIGHT_PHYSICSMODEL", "1")
    env["DRIVERFOAM_CARDIACFOAM_BACKEND"] = "lightweight"
    return env


def _driverfoam(*args: str, timeout: int = 300) -> subprocess.CompletedProcess:
    argv = [sys.executable, "-m", "openfoam_driver", *args]
    return subprocess.run(
        argv, env=_driver_env(), capture_output=True, text=True, timeout=timeout,
    )


def _plan_strict(tutorials_root: Path, entry: str, *, entry_kind: str | None = None) -> dict:
    args = ["plan", "--strict", "--entry", entry, "--tutorials-root", str(tutorials_root)]
    if entry_kind:
        args += ["--entry-kind", entry_kind]
    proc = _driverfoam(*args)
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        pytest.fail(
            f"driverFoam plan --strict produced no JSON for entry={entry!r}:\n"
            f"stdout={proc.stdout!r}\nstderr={proc.stderr[-2000:]}"
        )


def _plan_strict_tolerating_known_fallback_bug(tutorials_root: Path, entry: str) -> dict:
    """`_plan_strict`, tightened after the personalizedTemplates catalog
    fix (see module docstring's "FIXED BUG" section). This used to tolerate
    one specific false-positive diagnostic
    (ecgDomains.<name>.personalizedTemplates.ionicModelConfig.ionicModel
    spuriously "required" for a plain eikonalECG domain with no
    personalizedTemplates block); now that dict_entries_catalog.py no
    longer carries a typical_value on nBeats/duration/dt, dict-synthesis no
    longer ghost-populates them for domains that never configured
    personalizedTemplates, and that diagnostic no longer fires. This now
    requires a clean "ok" plan, same as `_plan_strict` -- kept as a
    separate name (rather than folded back into `_plan_strict`) so its
    callers document that they used to need this leniency."""
    plan = _plan_strict(tutorials_root, entry)
    assert plan.get("status") == "ok", (
        "plan --strict failed for the fallback (personalizedTemplates-"
        "stripped) case -- this should now plan cleanly after the "
        "personalizedTemplates catalog typical_value fix "
        f"(dict_entries_catalog.py): {plan.get('validation_diagnostics')!r}"
    )
    return plan


def _run_from_plan(plan: dict, run_doc_path: Path) -> dict:
    run_doc_path.write_text(json.dumps(plan["run_document"]))
    proc = _driverfoam("run", "--run-document", str(run_doc_path))
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        pytest.fail(
            f"driverFoam run --run-document produced no JSON:\n"
            f"stdout={proc.stdout!r}\nstderr={proc.stderr[-2000:]}"
        )


# --------------------------------------------------------------------------- #
# Case staging (copied temporary cases only -- never the tracked tutorial).
# --------------------------------------------------------------------------- #

_TUTORIAL_REL = "electrophysiologyProtocols/eikonalECGPersonalized"


def _stage_case(tmp_path_factory, *, tag: str) -> tuple[Path, Path]:
    root = tmp_path_factory.mktemp(f"eikonalecg_{tag}") / "tutorials"
    case_path = root / _TUTORIAL_REL
    case_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(monorepo_root / "tutorials" / _TUTORIAL_REL, case_path)
    return root, case_path


def _electro_properties(case_path: Path) -> Path:
    return case_path / "constant" / "electroProperties"


def _strip_personalized_templates_block(path: Path) -> None:
    """Remove the whole `personalizedTemplates { ... }` block via balanced-
    brace matching, leaving a plain eikonalECG domain (compiled-fallback
    path). Matches only the real key (`personalizedTemplates\\s*\\{`), not
    the word appearing in the file's own header comment."""
    text = path.read_text()
    match = re.search(r"personalizedTemplates\s*\n\s*\{", text)
    assert match, f"personalizedTemplates block not found in {path}"
    brace_start = text.index("{", match.end() - 1)
    depth = 0
    i = brace_start
    while True:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    end = i + 1
    line_start = text.rfind("\n", 0, match.start()) + 1
    path.write_text(text[:line_start] + text[end:].lstrip("\n"))


def _add_gkr_scale_override(path: Path, factor: float) -> None:
    """Fill in the tutorial's (committed, empty) `global { }`
    ionicConstantOverrides slot with a scale override on AC_GKr -- TWorld's
    rapid delayed-rectifier potassium conductance (src/ionicModels/
    TWorldBatched/TWorld_2025Batch.H:476, IKr current). See the module-level
    verification note in test_parameter_sensitivity_gkr_scaling_changes_ecg_
    trace for why this constant (spelled AC_GKr, not GKr -- confirmed via
    the solver's own "Unknown ionic constant" error listing) was chosen."""
    text = path.read_text()
    assert text.count("global { }") == 1, (
        f"expected exactly one empty 'global {{ }}' override slot in {path}"
    )
    path.write_text(
        text.replace("global { }", f"global {{ scale {{ AC_GKr {factor}; }} }}")
    )


@dataclass(frozen=True)
class CaseResult:
    case_path: Path
    plan: dict
    run: dict
    ecg: np.ndarray  # columns: time, E1, E2, E3, E4
    log_text: str


def _drive_case(
    tmp_path_factory, *, tag: str, mutate=None, tolerate_known_bug: bool = False,
) -> CaseResult:
    tutorials_root, case_path = _stage_case(tmp_path_factory, tag=tag)
    if mutate is not None:
        mutate(_electro_properties(case_path))

    if tolerate_known_bug:
        plan = _plan_strict_tolerating_known_fallback_bug(tutorials_root, _TUTORIAL_REL)
    else:
        plan = _plan_strict(tutorials_root, _TUTORIAL_REL)
        assert plan.get("status") == "ok", (
            f"plan --strict failed unexpectedly for {tag}: "
            f"{plan.get('validation_diagnostics')}"
        )

    run = _run_from_plan(plan, case_path.parent / "run_document.json")
    assert run.get("status") == "ok", f"run --run-document failed for {tag}: {run}"
    assert run.get("workflow_state", {}).get("status") == "completed", (
        f"workflow did not complete for {tag}: {run.get('workflow_state')}"
    )

    ecg_path = case_path / "postProcessing" / "eikonalECG.dat"
    assert ecg_path.is_file(), f"{ecg_path} was not written for {tag}"
    ecg = np.loadtxt(ecg_path, comments="#")

    log_path = case_path / "log.cardiacFoam"
    log_text = log_path.read_text() if log_path.is_file() else ""

    return CaseResult(case_path=case_path, plan=plan, run=run, ecg=ecg, log_text=log_text)


# --------------------------------------------------------------------------- #
# Module-scoped fixtures -- each real driverFOAM run happens exactly once
# and is shared across the checks that need it.
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def personalized_result(tmp_path_factory) -> CaseResult:
    """The eikonalECGPersonalized tutorial exactly as committed:
    personalizedTemplates present, ionicConstantOverrides.global empty."""
    return _drive_case(tmp_path_factory, tag="personalized")


@pytest.fixture(scope="module")
def fallback_result(tmp_path_factory) -> CaseResult:
    """Same tutorial with the personalizedTemplates block removed -> the
    compiled tissueTemplates.H fallback path. Previously had to tolerate
    one known catalog false-positive (module docstring's "FIXED BUG"
    section); now requires a clean "ok" plan like every other case here."""
    return _drive_case(
        tmp_path_factory,
        tag="fallback",
        mutate=_strip_personalized_templates_block,
        tolerate_known_bug=True,
    )


@pytest.fixture(scope="module")
def gkr_scaled_result(tmp_path_factory) -> CaseResult:
    """Same tutorial as `personalized_result`, but AC_GKr scaled to 30% in
    ionicConstantOverrides.global -- a real, verified-effective TWorld
    constant (see test_parameter_sensitivity_gkr_scaling_changes_ecg_trace)."""
    return _drive_case(
        tmp_path_factory,
        tag="gkr_scaled",
        mutate=lambda p: _add_gkr_scale_override(p, 0.3),
    )


# --------------------------------------------------------------------------- #
# Check 1 -- Fallback parity: the pre-existing manufacturedEikonalECG
# regression (no personalizedTemplates block) is unaffected.
# --------------------------------------------------------------------------- #

def test_fallback_parity_existing_manufactured_eikonal_ecg_regression_still_reproduces():
    """manufacturedSolutions/eikonalECG has no personalizedTemplates block,
    so its committed regression/eikonalECG.reference must still match
    byte-for-byte (within the reference's own declared tolerances) after
    this task's additions. Driven via the "generic" driver -- the case's
    own committed regression/regressionTest.sh -- which bypasses
    driverFOAM's Python dict-synthesis/catalog layer entirely, so it is
    unaffected by the known catalog bug documented in the module
    docstring. This is the strongest available "unchanged" check: it
    reruns the actual solver against the actual pre-existing reference."""
    from openfoam_driver.tests.regression_equivalence.dual_run import (
        solver_available, verify_reproduction,
    )
    from openfoam_driver.tests.regression_equivalence.registry import REGRESSION_CASES

    assert solver_available(), "cardiacFoam/WM_PROJECT_DIR unavailable despite requires_openfoam"

    case = next(c for c in REGRESSION_CASES if c.entry_name == "manufacturedEikonalECG")
    result = verify_reproduction(case, driver="generic")

    assert result.status == "reproduced", (
        f"manufacturedEikonalECG regression no longer reproduces via the "
        f"generic driver: {result.detail}"
    )


# --------------------------------------------------------------------------- #
# Check 2 -- Opt-in smoke test.
# --------------------------------------------------------------------------- #

_GENERATOR_SOURCE = (
    monorepo_root / "src/electroModels/ecgModels/eikonalECG/eikonalTemplateGenerator.C"
    if monorepo_root is not None else None
)


def test_opt_in_smoke_completes_writes_finite_ecg_and_generates_exactly_three_templates(
    personalized_result,
):
    # --- Static half: the generator's contract is fixed at exactly three
    # named anchors, each gated behind validateTemplate() -- which raises a
    # FatalError (aborting the whole run before any ECG sampling) if that
    # anchor's generated trace is invalid. This only reads the three anchor
    # names/count out of the real source; it does not reimplement or
    # second-guess validateTemplate()'s own logic.
    source = _GENERATOR_SOURCE.read_text()
    anchors = re.findall(r'validateTemplate\([^,]+,\s*"([^"]+)"\)', source)
    assert anchors == ["endocardium", "mid-myocardium", "epicardium"], (
        "eikonalTemplateGenerator.C's validated-anchor set changed from the "
        "three this smoke test's 'exactly three templates' claim depends on"
    )

    # --- Dynamic half: the run actually reached the post-generation ECG
    # sampling stage. Because generatePersonalizedTemplates() (which calls
    # validateTemplate() on all three anchors) runs synchronously, inside
    # solve(), strictly before the "wrote sampled ECG" log line, reaching
    # that line is direct evidence all three anchors were generated and
    # passed validation -- not merely that the case didn't crash for some
    # unrelated reason.
    assert "eikonalECG: wrote sampled ECG" in personalized_result.log_text

    ecg = personalized_result.ecg
    assert np.isfinite(ecg).all(), "eikonalECG.dat contains non-finite values"
    # sampling.start=0, sampling.end=0.6, sampling.deltaT=1e-4 -> 6001 rows,
    # 5 columns (time + 4 electrodes E1..E4).
    assert ecg.shape == (6001, 5)
    assert ecg[0, 0] == pytest.approx(0.0)
    assert ecg[-1, 0] == pytest.approx(0.6)


# --------------------------------------------------------------------------- #
# Check 3 -- Baseline parity against the compiled fallback.
#
# Investigation (see task-3b-report.md for the full numbers): running both
# paths on identical copies of the tutorial and diffing postProcessing/
# eikonalECG.dat shows the two traces are effectively IDENTICAL for
# t >= 0.01s (correlation 0.9999, RMSE ~3e-6, peak amplitudes agreeing to 5
# significant figures) -- this is the activation-time-field-driven
# repolarization/T-wave plateau, which both paths reconstruct from
# fundamentally the same eikonal solve and the same ionicHeterogeneity
# blend weights. All of the discrepancy is concentrated in the first ~1ms
# (the QRS-like upstroke transient), where the personalized path's peak is
# roughly 4.3x the compiled fallback's. Since the pseudo-ECG source term is
# gradVm ~ dV/dt (not V itself), this is exactly the signature of the two
# templates' upstrokes having different original sampling/pacing history --
# precisely the "not recoverable, expected" divergence the task brief
# anticipated (no tissueTemplates.H generation script exists in this repo;
# git log --follow on it returns a single squashed commit) -- not a wrong-
# shape/wrong-sign/order-of-magnitude structural bug: both peaks are
# positive, occur within a fraction of a millisecond of each other, and the
# ratio (4.3x) is nowhere near the ~1e2-1e3x a genuine unit/scaling bug in
# the mV/V conversion would produce.
#
# The tolerance below is therefore two-part, set from what was actually
# observed (not guessed): a tight bound for t >= 0.01s (>15x margin over the
# observed ~3e-6 RMSE there), and a loose-but-bounded envelope for the
# upstroke transient that still catches a real bug (wrong sign, a zero
# trace, or an order-of-magnitude-plus blowup) without failing on the
# expected, documented amplitude difference.
# --------------------------------------------------------------------------- #

_LATE_PHASE_START_S = 0.01
_LATE_PHASE_ATOL = 5e-5   # V; ~15x the observed ~3e-6 RMSE there
_LATE_PHASE_RTOL = 0.05
_EARLY_PEAK_RATIO_BOUNDS = (0.05, 20.0)  # observed ~4.3x; ample margin either side


def test_baseline_parity_against_compiled_fallback(personalized_result, fallback_result):
    pers, fb = personalized_result.ecg, fallback_result.ecg
    assert pers.shape == fb.shape
    t = pers[:, 0]
    assert np.allclose(t, fb[:, 0]), "the two runs sampled different time grids"

    late = t >= _LATE_PHASE_START_S
    early = ~late
    assert early.any() and late.any()

    for col, name in enumerate(("E1", "E2", "E3", "E4"), start=1):
        p_late, f_late = pers[late, col], fb[late, col]
        assert np.allclose(p_late, f_late, atol=_LATE_PHASE_ATOL, rtol=_LATE_PHASE_RTOL), (
            f"{name}: post-upstroke (t>={_LATE_PHASE_START_S}s) traces diverge "
            "beyond the observed-noise-floor tolerance -- this is the region "
            "that should agree almost exactly between the two paths"
        )

        p_early, f_early = pers[early, col], fb[early, col]
        assert np.isfinite(p_early).all() and np.isfinite(f_early).all()
        p_peak_idx = np.argmax(np.abs(p_early))
        f_peak_idx = np.argmax(np.abs(f_early))
        p_peak, f_peak = p_early[p_peak_idx], f_early[f_peak_idx]
        assert p_peak != 0 and f_peak != 0, f"{name}: upstroke transient is exactly zero"
        assert np.sign(p_peak) == np.sign(f_peak), (
            f"{name}: personalized and fallback upstroke peaks have opposite sign "
            f"(pers={p_peak}, fallback={f_peak}) -- looks like a real bug, not "
            "expected pacing-history drift"
        )
        ratio = abs(p_peak) / abs(f_peak)
        lo, hi = _EARLY_PEAK_RATIO_BOUNDS
        assert lo <= ratio <= hi, (
            f"{name}: upstroke peak ratio {ratio:.3g} outside the envelope "
            f"[{lo}, {hi}] set from the observed ~4.3x difference -- this is "
            "large enough to suggest a real bug (unit/scaling error), not "
            "just unrecoverable pacing-history drift"
        )


# --------------------------------------------------------------------------- #
# Check 4 -- Parameter sensitivity (AC_GKr).
#
# Verification that AC_GKr is real and effective for TWorld (done BEFORE
# building this check, per the brief): an isolated single-cell run
# (tutorials/electrophysiologyProtocols/singleCell, ionicModel TWorld,
# copied and driven the same way as this file's cases) with
# ionicConstantOverrides.global.scale.AC_GKr=0.3 shifted APD90 from 269.0ms
# to 406.0ms (+137ms, using the tutorial's own calc_apd.py APD90 formula --
# threshold crossing at -10mV for upstroke, 90%-repolarization crossing for
# offset) relative to an override-free baseline run of the same case. This
# is a large, physiologically sane effect (reduced IKr prolongs
# repolarization) and confirms AC_GKr is a real, effective TWorld constant
# reachable through ionicConstantOverrides -- not an assumption. (Also
# confirms the correct spelling is "AC_GKr", not "GKr": the solver's own
# "Unknown ionic constant" FatalError lists AC_GKr/AC_GKr_b among TWorld's
# 274 real constant names.)
#
# The RMSE threshold below is set from what was actually observed on the
# eikonalECG trace itself: personalized-baseline vs. AC_GKr*0.3 gives a
# full-trace RMSE of ~2.2e-4 V, about 70x the baseline-parity check's
# late-phase noise floor (~3e-6 V) and far above the ~0 floor of rerunning
# the identical config twice (confirmed deterministic to 9 significant
# figures). The threshold (5e-5) sits with a comfortable margin on both
# sides of that separation.
# --------------------------------------------------------------------------- #

_SENSITIVITY_RMSE_THRESHOLD = 5e-5  # V; observed effect ~2.2e-4, noise floor ~3e-6


def test_parameter_sensitivity_gkr_scaling_changes_ecg_trace(personalized_result, gkr_scaled_result):
    pers, scaled = personalized_result.ecg, gkr_scaled_result.ecg
    assert pers.shape == scaled.shape
    assert np.allclose(pers[:, 0], scaled[:, 0])

    diff = scaled[:, 1:] - pers[:, 1:]
    rmse = float(np.sqrt(np.mean(diff**2)))
    assert rmse > _SENSITIVITY_RMSE_THRESHOLD, (
        f"scaling AC_GKr to 0.3x produced too small an ECG change "
        f"(rmse={rmse:.3g}, threshold={_SENSITIVITY_RMSE_THRESHOLD:.3g}) -- "
        "either the override isn't reaching the generator, or AC_GKr no "
        "longer measurably affects TWorld's repolarization here"
    )
    assert np.isfinite(scaled).all()


# --------------------------------------------------------------------------- #
# Check 5 -- Unit guard.
#
# Deterministic, ionic-model-amplitude-independent check of the shared
# piecewise-linear derivative primitive (tissueTemplates.H::
# evaluateTemplateDerivative, used by all three of eikonalECG.C's
# reconstructGradVm branches -- manufactured, personalized, compiled
# fallback), the mV/s -> V/s conversion, and the chain-rule minus sign
# (eikonalECG.C: "Vm(x,t) = U(t-tau(x)) => gradVm = -dU/ds(t-tau) * gradTau").
#
# A tiny throwaway program (same technique as src/electroModels/tests/
# test_eikonal_template_generator.py -- tissueTemplates.H is header-only
# via scalar.H, so no linking against libelectroModels is actually needed
# here) evaluates a synthetic two-point template (0mV at t=0s, 1mV at
# t=0.001s -> an exact secant slope of 1000 mV/s) at a point strictly
# inside the segment, then applies eikonalECG.C's own conversion+chain-rule
# lines verbatim. A separate static check confirms those exact lines are
# still present, unchanged, in the real solver source.
# --------------------------------------------------------------------------- #

_EIKONAL_ECG_SOURCE = (
    monorepo_root / "src/electroModels/ecgModels/eikonalECG/eikonalECG.C"
    if monorepo_root is not None else None
)
_TISSUE_TEMPLATES_INCLUDE_DIR = (
    monorepo_root / "src/electroModels/ecgModels/eikonalECG"
    if monorepo_root is not None else None
)

_UNIT_GUARD_CPP = '''\
#include "tissueTemplates.H"
#include <cstdio>

using namespace Foam;
using namespace Foam::eikonalECG_templates;

int main()
{
    // Two-point synthetic template: exact secant slope of 1000 mV/s,
    // sampled strictly inside the segment (not at either endpoint, which
    // would hit evaluateTemplateDerivative's clamped-zero branches).
    static const scalar times[]  = {0.0, 0.001};
    static const scalar values[] = {0.0, 1.0};

    const scalar rawDUds = evaluateTemplateDerivative(0.0005, times, values, 2);

    // Verbatim reproduction of eikonalECG.C's reconstructGradVm conversion
    // + chain rule ("Dynamic values are mV, so convert dVm/dt to V/s
    // before the chain rule." / "Vm(x,t) = U(t-tau(x)) =>
    // gradVm = -dU/ds(t-tau) * gradTau"), with a unit activation-time
    // gradient standing in for gradTauValues[cellI].
    const scalar dUds = rawDUds*1e-3;
    const scalar unitGradTau = 1.0;
    const scalar gradVmComponent = -dUds * unitGradTau;

    std::printf("%.10f %.10f %.10f\\n", rawDUds, dUds, gradVmComponent);
    return 0;
}
'''


def _build_unit_guard_program(tmp_path: Path) -> Path:
    compiler = _compiler()
    assert compiler is not None  # gated by requires_openfoam's C++ toolchain expectation

    cpp_file = tmp_path / "unit_guard.cpp"
    exe_file = tmp_path / "unit_guard"
    cpp_file.write_text(_UNIT_GUARD_CPP)

    wm_project_dir = Path(os.environ["WM_PROJECT_DIR"])
    compile_cmd = [
        compiler,
        "-std=c++17", "-m64", "-pthread", "-ftrapping-math",
        "-DOPENFOAM=2412", "-DWM_DP", "-DWM_LABEL_SIZE=32",
        "-O3", "-DNoRepository", "-ftemplate-depth-100",
        "-DOPENFOAM_COM", "-DOPENFOAM_NOT_EXTEND",
        "-I", str(_TISSUE_TEMPLATES_INCLUDE_DIR),
        "-I", str(wm_project_dir / "src/OpenFOAM/lnInclude"),
        "-I", str(wm_project_dir / "src/OSspecific/POSIX/lnInclude"),
        "-fPIC",
        str(cpp_file), "-o", str(exe_file),
    ]
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"unit-guard program failed to compile:\n{result.stderr}"
    )
    return exe_file


def test_unit_guard_two_point_derivative_and_chain_rule_sign(tmp_path):
    exe = _build_unit_guard_program(tmp_path)
    result = subprocess.run([str(exe)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr

    raw_dUds, dUds, grad_component = (float(x) for x in result.stdout.split())

    # 1000 mV/s from the shared interpolation primitive.
    assert raw_dUds == pytest.approx(1000.0, abs=1e-9)
    # 1 V/s after the documented *1e-3 conversion -- BEFORE the chain-rule
    # minus sign, per the brief's exact wording.
    assert dUds == pytest.approx(1.0, abs=1e-12)
    # Sign-flipped by the chain rule against a unit activation-time gradient.
    assert grad_component == pytest.approx(-1.0, abs=1e-12)

    # Static corroboration: the literal conversion factor and chain-rule
    # sign this program reproduces verbatim still appear, unchanged, in the
    # real solver source -- so a change to either literal is caught here
    # even though this program only mirrors, rather than calls, that code.
    source = _EIKONAL_ECG_SOURCE.read_text()
    assert "dUds = rawDUds*1e-3;" in source
    assert "gradVmValues[cellI] = -dUds * gradTauValues[cellI];" in source
