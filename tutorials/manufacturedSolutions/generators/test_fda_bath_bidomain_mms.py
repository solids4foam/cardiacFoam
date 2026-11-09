import importlib.util
import os
import re

import numpy as np
import pytest
import sympy as sp

import fda_bath_bidomain_mms as fbb
import fda_bidomain_mms as fbm
import fda_monodomain_mms as fmm
import slab_leadfield_mms as slm

DIMS = (1, 2, 3)
x, y, z, t = fbb.x, fbb.y, fbb.z, fbb.t
k, alpha, se = fbb.k, fbb.alpha, fbb.se
CASE = os.path.join(fbb.HERE, "..", "bathBidomain")
IONIC = os.path.join(fbb.HERE, "..", "..", "..", "src", "ionicModels", "verificationModels",
                     "bathBidomainFDAManufactured", "bathBidomainFDAManufactured_2014.H")


def _zero(e):
    return sp.simplify(e) == 0


def test_self_checks_and_cli():
    assert fbb.main(["--self-check"]) == 0


def test_constants():
    assert (fbb.CM, fbb.CHI) == (2, 3)
    assert fbb.GI_DIAG == fbm.GI_DIAG and fbb.GE_DIAG == fbm.GE_DIAG
    assert fbb.SIGMA_B == sp.Rational("0.0230827345") == fbb.GE_DIAG[0]/2
    assert fbb.ALPHA == sp.Rational(1, 100)
    assert fbb.K == 1/sp.sqrt(2)
    assert fbb.BETA_SCALE == sp.Rational(-11, 10)
    assert set(fbb.fields(3)) == {"Vm", "phiE", "phiI", "u1", "u2", "u3", "Iion", "F", "G",
                                  "C", "sigmaB"}


@pytest.mark.parametrize("d", DIMS)
def test_ionic_rates_hold(d):
    f = fbb.fields(d)
    rates = fmm.ionic_rates(f["Vm"], f["u1"], f["u2"], f["u3"])
    for name in ("u1", "u2", "u3"):
        assert _zero(sp.diff(f[name], t) - rates[name]), name


@pytest.mark.parametrize("d", DIMS)
def test_heart_residuals_are_rounding_level(d):
    """V_m equation with beta = -1.1 (1 - k) and the elliptic equation, at the
    tutorial's k, relative to Gi_xx pi^2: rounding of the nine-decimal input."""
    rel_p, rel_e = fbb.relative_residuals(d)
    assert rel_p <= fbb.RESIDUAL_RTOL and rel_e <= fbb.RESIDUAL_RTOL


def test_interfaces_and_boundaries():
    branches = fbb.phie_branches()
    sb = fbb.sigma_b_of(se)
    left, heart, right = (b[0].subs(fbb.SB, sb) for b in branches)
    phiI = fbb.fields(1)["Vm"] + heart
    for wall, a, b in ((0, left, heart), (1, heart, right)):
        assert _zero((a - b).subs(x, wall)), "phi_e continuity"
        bath = left if wall == 0 else right
        assert _zero((se*sp.diff(heart, x) - sb*sp.diff(bath, x)).subs(x, wall)), "current"
        assert _zero(sp.diff(phiI, x).subs(x, wall)), "sealed intracellular"
    assert _zero(left.subs(x, -1)), "ground on x = -1"
    assert _zero(sb*sp.diff(right, x) - alpha), "surface current on x = 2"
    for bath in (left, right):
        assert _zero(sp.diff(bath, x, 2))


def test_alpha_zero_is_the_slab():
    phiE = fbb.fields(1)["phiE"].subs(alpha, 0)
    slab = slm.fields()["phiSlab"].subs(slm.k, k)
    for xv in (-0.7, 0, 0.3, 1, 1.6):
        assert _zero((phiE - slab).subs(x, xv))


def test_emitted_files_up_to_date():
    assert fbb.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    fbb.write_outputs(str(tmp_path))
    for rel, text in fbb.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference():
    path = os.path.join(fbb.HERE, "reference", "fda_bath_bidomain.py")
    spec = importlib.util.spec_from_file_location("fda_bath_bidomain_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("d", DIMS)
def test_emitted_numpy_matches_sympy(d):
    ref = _reference()
    rng = np.random.default_rng(7)
    pts = np.column_stack([rng.uniform(-1, 2, 8), rng.uniform(0, 1, (8, 2)),
                           rng.uniform(0, 0.1, 8)])
    par = (float(fbb.K), float(fbb.ALPHA), float(fbb.GE_DIAG[0]))
    for name in ("Vm", "phiE", "phiI", "u1", "u2", "u3"):
        expr = fbb.fields(d)[name].subs(fbb.SB, fbb.sigma_b_of(se))
        fn = sp.lambdify((x, y, z, t, k, alpha, se), expr, "numpy")
        want = np.broadcast_to(fn(*pts.T, *par), (8,))
        got = getattr(ref, "%s_d%d" % (name, d))(*pts.T, *par)
        np.testing.assert_allclose(got, want, rtol=1e-13, atol=1e-15)


@pytest.mark.skipif(not os.path.isdir(CASE), reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_matches():
    props = open(os.path.join(CASE, "constant", "electroProperties")).read()
    tensor = r"\s*\[[^\]]*\]\s*\(([^)]*)\)"
    gi = re.search(r"conductivityIntracellular" + tensor, props).group(1).split()
    ge = re.search(r"conductivityExtracellular" + tensor, props).group(1).split()
    assert [sp.Rational(v) for v in gi] == [fbb.GI_DIAG[0], 0, 0, fbb.GI_DIAG[1], 0, fbb.GI_DIAG[2]]
    assert [sp.Rational(v) for v in ge] == [fbb.GE_DIAG[0], 0, 0, fbb.GE_DIAG[1], 0, fbb.GE_DIAG[2]]
    assert float(re.search(r"\bk\s+(\S+);", props).group(1)) == float(fbb.K)
    assert sp.Rational(re.search(r"\balpha\s+(\S+);", props).group(1)) == fbb.ALPHA
    assert re.search(r"\bchi\s*\[[^\]]*\]\s*(\S+);", props).group(1) == "3"
    assert re.search(r"\bcm\s*\[[^\]]*\]\s*(\S+);", props).group(1) == "2"
    cond = open(os.path.join(CASE, "system", "setTorsoOrganConductivityFieldDict")).read()
    assert sp.Rational(re.search(r"\bbath\s+(\S+);", cond).group(1)) == fbb.SIGMA_B
    ionic = open(IONIC).read()
    assert re.search(r"betaScale\s*=\s*(\S+);", ionic).group(1) == "-1.1"
    assert "CONSTANTS[BathBidomainCm] = 2.0;" in ionic
    assert "CONSTANTS[BathBidomainChi] = 3.0;" in ionic


# ---------------------------------------------------- C++ header regions ----

@pytest.mark.parametrize("name", sorted(fbb.HEADER_PATHS))
def test_header_up_to_date(name):
    assert fbb.header_matches(fbb.HEADER_PATHS[name])


@pytest.mark.parametrize("name", sorted(fbb.HEADER_PATHS))
def test_header_check_cli(name):
    assert fbb.main(["--check", fbb.HEADER_PATHS[name]]) == 0


def test_header_write_touches_only_regions(tmp_path):
    import mms_codegen as codegen
    for name, rel in fbb.HEADER_PATHS.items():
        text = open(codegen.resolve(fbb.__file__, rel)).read()
        stale = codegen.splice_regions(text, {n: ["// stale"] for n in fbb.HEADERS[name](text)})
        assert stale != text
        dst = tmp_path/name
        dst.write_text(stale)
        assert fbb.main(["--check", str(dst)]) == 1
        assert fbb.main(["--write", str(dst)]) == 0
        assert dst.read_text() == text


def test_header_without_regions_is_refused(tmp_path):
    for name in fbb.HEADER_PATHS:
        dst = tmp_path/name
        dst.write_text("namespace Foam {}\n")
        assert fbb.main(["--check", str(dst)]) == 1
        assert fbb.main(["--write", str(dst)]) == 1
        assert dst.read_text() == "namespace Foam {}\n"


@pytest.mark.parametrize("d", DIMS)
def test_parity_bath(d):
    """The compiled header against reference/fda_bath_bidomain.py and the
    sympy gradients at the tutorial's k, alpha and se."""
    try:
        from parity.test_parity import compare, parity_available, sample_points
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()
    p = sample_points()
    p[:, 0] = 3.0*p[:, 0] - 1.0
    # not x = -1: G_1 = 1 + x vanishes there and u_2 is infinite (the ionic
    # states live in the heart; phi_e(-1) = 0 is a self-check)
    p = np.vstack([p, [[v, 0.37, 0.58, 0.03] for v in (-0.999, 0.0, 1.0, 2.0)]])
    par = (float(fbb.K), float(fbb.ALPHA), float(fbb.GE_DIAG[0]))
    sub = {fbb.SB: fbb.sigma_b_of(se)}
    heart = fbb.phie_branches()[1][0].subs(sub)
    gV = sp.lambdify((x, t, k, alpha, se), sp.diff(fbb.fields(1)["Vm"], x), "numpy")
    gE = sp.lambdify((x, t, k, alpha, se), sp.diff(heart, x), "numpy")
    gI = sp.lambdify((x, t, k, alpha, se), sp.diff(fbb.fields(1)["Vm"] + heart, x), "numpy")
    gB = float(fbb.ALPHA/fbb.SIGMA_B)

    def numpy_fn(q):
        a = [q[:, i] for i in range(4)]
        n = q.shape[0]
        cols = [getattr(ref, "%s_d%d" % (m, d))(*a, *par)
                for m in ("Vm", "u1", "u2", "u3", "phiE", "phiI")]
        cols += [getattr(ref, "%s_d%d" % (m, d))(*a, *par) for m in ("Vm", "phiE", "phiI")]
        cols += [np.broadcast_to(gV(q[:, 0], q[:, 3], *par), (n,)),
                 np.broadcast_to(gE(q[:, 0], q[:, 3], *par), (n,)),
                 np.broadcast_to(gI(q[:, 0], q[:, 3], *par), (n,)), np.full(n, gB)]
        return np.column_stack(cols)

    compare("bathBidomain", d, numpy_fn, p, rtol=1e-12,
            k="%.17g" % par[0], alpha="%.17g" % par[1], se="%.17g" % par[2])


METRICS = os.path.join(fbb.HERE, "..", "..", "..", "applications", "utilities",
                       "bathBidomainInterfaceMetrics", "bathBidomainInterfaceMetrics.C")


@pytest.mark.skipif(not os.path.exists(METRICS), reason="cardiacFoam-only")
def test_interface_metrics_uses_the_generated_potential():
    """bathBidomainInterfaceMetrics evaluates the exact phi_e through the
    generated header, not a copy of its branches."""
    text = open(METRICS).read()
    assert '#include "bathBidomainVerification/manufacturedFDABathBidomainReference.H"' in text
    assert "manufacturedFDABathPhiEAt(" in text
    assert not re.search(r"\bscalar\s+phiEExact\s*\(", text), "hand-written phiEExact"
    assert "sqrt1t" not in text and "0.5*sigmaE" not in text
