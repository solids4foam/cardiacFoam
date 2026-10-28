import importlib.util
import os

import numpy as np
import pytest
import sympy as sp

import fda_monodomain_mms as fda
import graph_monodomain_mms as gr

x, t = gr.x, gr.t


def test_graph_residual_is_rounding_level():
    # No source is applied on the graph: the residual of the equation is
    # (g pi^2 + beta_1) V_1D, below 1e-8 |beta_1| |V_1D|.
    f = gr.fields()
    assert sp.simplify(gr.graph_residual() - f["residual_1D"]) == 0
    assert "S1D" not in f
    k1 = gr.graph_conductivity()*sp.pi**2
    assert sp.simplify(f["residual_1D"] - (k1 + sp.Rational(str(fda.BETA[1])))*f["V1D"]) == 0
    assert abs(float(k1) + fda.BETA[1]) <= 1e-8*abs(fda.BETA[1])


def test_graph_constants():
    assert gr.G_GRAPH == 0.111453302 == float(fda.SIGMA_DIAG[0])
    assert gr.ROOT == (sp.Rational(1, 2), sp.Rational(1, 6), sp.Rational(1, 3))
    assert gr.TERMINALS[0] == (0, sp.Rational(1, 6), sp.Rational(1, 3))
    assert gr.TERMINALS[1] == (1, sp.Rational(1, 6), sp.Rational(1, 3))


def test_graph_field_is_fda_d1():
    f = gr.fields()
    assert sp.simplify(f["V1D"] - sp.sqrt(1 + t)*sp.cos(sp.pi*x)) == 0
    for n in ("u1", "u2", "u3"):
        assert sp.simplify(f[n] - fda.fields(1)[n]) == 0


def test_terminal_values():
    Vm = gr.tissue_vm()
    V1 = gr.fields()["V1D"]
    X, Y, Z = fda.x, fda.y, fda.z
    for tx, ty, tz in gr.TERMINALS:
        at = {X: tx, Y: ty, Z: tz}
        assert sp.simplify(Vm.subs(at) + V1.subs(x, tx)/2) == 0
        assert sp.simplify(sp.diff(Vm, X).subs(at)) == 0


def test_self_checks_and_cli():
    assert gr.main(["--check"]) == 0


def test_emitted_files_up_to_date():
    assert gr.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    gr.write_outputs(str(tmp_path))
    for rel, text in gr.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference():
    path = os.path.join(gr.HERE, "reference", "graph_monodomain.py")
    spec = importlib.util.spec_from_file_location("graph_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_emitted_numpy_matches_sympy():
    ref = _reference()
    pts = np.random.default_rng(6).uniform(0, 1, (6, 2))
    for name, expr in gr.fields().items():
        fn = sp.lambdify((x, t), expr, "numpy")
        want = np.broadcast_to(fn(*pts.T), (6,))
        np.testing.assert_allclose(getattr(ref, name)(*pts.T), want, rtol=1e-13, atol=1e-15)


def test_parity_graph():
    try:
        from parity.test_parity import compare, parity_available, sample_points
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    ref = _reference()

    def numpy_fn(p):
        a = (p[:, 0], p[:, 3])
        return np.column_stack([getattr(ref, n)(*a) for n in ("V1D", "u1", "u2", "u3")])

    compare("graph", 1, numpy_fn, sample_points(), rtol=1e-12)
