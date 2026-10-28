"""Parity harness: the current hand-written C++ reference functions, evaluated
by the `mmsParity` OpenFOAM app, against a NumPy reference.

Only the cardiacFoam checkout can build `mmsParity`; a release folder has no
OpenFOAM, so generator tests import `compare` lazily and skip when
`parity_available()` reports why the harness cannot run.  Run directly
(`pytest parity/test_parity.py`) the smoke tests FAIL, not skip, without
OpenFOAM.

Protocol of `mmsParity <case> <d> [key=value ...]`
--------------------------------------------------
stdin   one point per line: `x y z t`.
stdout  one line per point: space-separated values with 17 significant digits.
        For `pseudoECG*` and `eikonalECG` the point `x y z` is the electrode
        position (the integral does not depend on a field point) and `t` the
        time.
keys    beta=<v>              ionic constant Beta (default per d: -1.1, -5.9, -8.6)
        order=<n>             Gauss-Legendre order for the ECG cases
                              (default 96 pseudoECG*, 48 eikonalECG)
        conductivity=rotated|axis   tensor of `eikonal`/`eikonalECG`
                                    (default rotated)
        c0=<v>                eikonal speed constant (default 1)
        template=singleEvent|sine   eikonalECG template (default singleEvent)
        duration=<v>, rise=<v>      singleEvent parameters (default 64, 4)

case              d        columns
fda               1,2,3    V u1 u2 u3 Iion
rotated           3        V u1 u2 u3 pdeSource
pseudoECG         1,2,3    spatialReference reference(t)   (G = diag FDA)
pseudoECGRotated  3        spatialReference reference(t)   (G = full tensor)
eikonal           1,2,3    tau gradTau_x gradTau_y gradTau_z sourceTerm
eikonalECG        1,2,3    reference(t)
graph             1,2,3    V u1 u2 u3 (graph node at x y z)
coupled           3        V1D_0 Vm3D_0 I_0 V1D_1 Vm3D_1 I_1 (terminals of
                           the case, R_pvj = 1; x y z ignored, t used)
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
APP = HERE / "mmsParity"
SOURCES = [HERE / "mmsParity.C", HERE / "Make" / "files", HERE / "Make" / "options"]
DEFAULT_BASHRC = Path("/Volumes/OpenFOAM-v2412/etc/bashrc")

# case -> (dimensions, number of output columns)
CASES = {
    "fda": ((1, 2, 3), 5),
    "rotated": ((3,), 5),
    "pseudoECG": ((1, 2, 3), 2),
    "pseudoECGRotated": ((3,), 2),
    "eikonal": ((1, 2, 3), 5),
    "eikonalECG": ((1, 2, 3), 1),
    "graph": ((1, 2, 3), 4),
    "coupled": ((3,), 6),
}


def _bashrc():
    return Path(os.environ.get("CARDIAC_OPENFOAM_BASHRC", DEFAULT_BASHRC))


def parity_available():
    """(True, "") when mmsParity can be built and run, else (False, reason)."""
    if not (HERE / "mmsParity.C").exists():
        return False, "parity harness absent (cardiacFoam-only)"
    if shutil.which("bash") is None:
        return False, "bash not found"
    if "WM_PROJECT" not in os.environ and not _bashrc().exists():
        return False, (
            f"OpenFOAM environment unavailable: WM_PROJECT is unset and "
            f"{_bashrc()} does not exist (set CARDIAC_OPENFOAM_BASHRC)"
        )
    return True, ""


def _bash(command):
    return subprocess.run(
        ["bash", "-c", f'source "{_bashrc()}" >/dev/null 2>&1; {command}'],
        capture_output=True, text=True,
    )


def build_app():
    """wmake mmsParity into parity/ when missing or older than its sources."""
    if APP.exists() and all(
        APP.stat().st_mtime >= s.stat().st_mtime for s in SOURCES
    ):
        return APP
    done = _bash(f'cd "{HERE}" && wmake')
    if done.returncode != 0 or not APP.exists():
        raise RuntimeError(
            f"wmake of mmsParity failed\n{done.stdout[-3000:]}\n{done.stderr[-3000:]}"
        )
    return APP


def sample_points(n=20, seed=20260101):
    """n seeded points: x, y, z uniform in [0, 1], t uniform in [0, 0.1]."""
    rng = np.random.default_rng(seed)
    return np.column_stack([rng.uniform(0, 1, (n, 3)), rng.uniform(0, 0.1, n)])


def run_app(case, d, points, **keys):
    """Evaluate the C++ reference at `points` (N, 4); returns (N, ncols)."""
    build_app()
    args = [str(APP), case, str(d)] + [f"{k}={v}" for k, v in keys.items()]
    stdin = "".join(
        " ".join(f"{v:.17g}" for v in row) + "\n" for row in np.atleast_2d(points)
    )
    done = _bash(
        " ".join(f"'{a}'" for a in args) + " <<'EOF_POINTS'\n" + stdin + "EOF_POINTS"
    )
    if done.returncode != 0:
        raise RuntimeError(f"mmsParity {case} {d} failed\n{done.stderr[-3000:]}")
    return np.array(
        [[float(v) for v in line.split()] for line in done.stdout.strip().splitlines()]
    )


def compare(case, d, numpy_fn, points, rtol=1e-12, **keys):
    """Assert numpy_fn(points) equals the C++ values.

    numpy_fn maps the (N, 4) array to (N,) or (N, ncols).  A value passes when
    |a - b| <= rtol |b|, or |a - b| <= rtol where |b| < 1e-8.
    """
    cpp = run_app(case, d, points, **keys)
    got = np.asarray(numpy_fn(points), dtype=float).reshape(cpp.shape)
    err = np.abs(got - cpp)
    tol = np.where(np.abs(cpp) < 1e-8, rtol, rtol * np.abs(cpp))
    bad = err > tol
    assert not bad.any(), (
        f"{case} d={d}: {bad.sum()} of {bad.size} values differ; worst "
        f"|diff|={err.max():.3e} at {np.unravel_index(err.argmax(), err.shape)}"
    )


@pytest.fixture(scope="module")
def app():
    ok, why = parity_available()
    if not ok:
        pytest.fail(why)
    return build_app()


CASE_D = [(c, d) for c, (ds, _) in CASES.items() for d in ds]


@pytest.mark.parametrize("case,d", CASE_D)
def test_smoke_finite_expected_count(app, case, d):
    out = run_app(case, d, sample_points())
    assert out.shape == (20, CASES[case][1])
    assert np.isfinite(out).all()


def test_smoke_single_point(app):
    out = run_app("fda", 3, np.array([[0.3, 0.2, 0.1, 0.05]]))
    assert out.shape == (1, 5) and np.isfinite(out).all()


def test_fda_closed_form_spot_value(app):
    # V = sqrt(1 + t) cos(pi x) at d = 1: pins the column order.
    out = run_app("fda", 1, np.array([[0.25, 0.0, 0.0, 0.5]]))
    assert out[0, 0] == pytest.approx(np.sqrt(1.5) * np.cos(np.pi * 0.25), rel=1e-14)


def test_electrode_on_quadrature_point_gives_zero_not_nan(app):
    # Order 1 has the single node x = 1/2: an electrode there sits on it.
    out = run_app("pseudoECG", 1, np.array([[0.5, 0.0, 0.0, 0.0]]), order=1)
    assert out[0, 0] == 0.0 and out[0, 1] == 0.0


def test_eikonal_ecg_template_key(app):
    pts = sample_points(3)
    a = run_app("eikonalECG", 3, pts, template="sine", order=8)
    b = run_app("eikonalECG", 3, pts, template="singleEvent", order=8)
    assert np.isfinite(a).all() and np.isfinite(b).all()
    assert not np.allclose(a, b)


def test_compare_detects_mismatch(app):
    with pytest.raises(AssertionError):
        compare("graph", 1, lambda p: np.zeros((len(p), 4)), sample_points(5))


def test_coupled_terminal_current_matches_continuum_value(app):
    # README of monodomain1D3D: with R_pvj = 1 the terminal currents are
    # +-(3/2) sqrt(1 + t); terminal 0 carries +, terminal 1 carries -.
    # The graph file holds y and z to 12 digits, hence rel 1e-10.
    t = 0.07
    out = run_app("coupled", 3, np.array([[0.0, 0.0, 0.0, t]]))
    s = np.sqrt(1 + t)
    assert out[0, 2] == pytest.approx(1.5 * s, rel=1e-10)
    assert out[0, 5] == pytest.approx(-1.5 * s, rel=1e-10)
