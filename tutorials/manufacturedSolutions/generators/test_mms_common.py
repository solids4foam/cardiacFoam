"""Tests for the shared emitters and quadratures in mms_common."""

import math
import sys

import numpy as np
import sympy as sp

from mms_common.latex_emit import PAPER_SYMBOLS, emit_latex_fragment, paper_symbol_names
from mms_common.numpy_emit import emit_numpy_module
from mms_common.quadrature import gauss_legendre_box, half_ball_kernel_integral


def test_numpy_module_imports_without_sympy(monkeypatch):
    x, y = sp.symbols("x y")
    src = emit_numpy_module("demo", {"f": ([x, y], sp.sin(x) * y)}, {"a": 1.5}, "Demo.")
    assert "sympy" not in src
    monkeypatch.setitem(sys.modules, "sympy", None)
    ns = {}
    exec(compile(src, "demo.py", "exec"), ns)
    assert ns["PARAMS"] == {"a": 1.5}
    assert abs(ns["f"](0.3, 2.0) - 2 * math.sin(0.3)) < 1e-15
    out = ns["f"](np.array([0.3, 0.4]), np.array([2.0, 1.0]))
    assert out.shape == (2,)


def test_numpy_module_matrix_and_constant_broadcast():
    x = sp.Symbol("x")
    src = emit_numpy_module(
        "demo", {"v": ([x], sp.Matrix([x**2, 1])), "c": ([x], sp.Integer(3))}, {}, "Demo.")
    ns = {}
    exec(compile(src, "demo.py", "exec"), ns)
    xs = np.array([1.0, 2.0, 3.0])
    v = ns["v"](xs)
    assert v.shape == (2, 3)
    assert np.allclose(v[0], xs**2) and np.allclose(v[1], 1.0)
    assert ns["c"](xs).shape == (3,)


def test_latex_uses_paper_symbols():
    expr = sp.Symbol("Vm") * sp.Symbol("chi")
    tex = emit_latex_fragment("a", expr, paper_symbol_names(expr.free_symbols))
    assert "V_m" in tex and r"\chi" in tex
    assert "Vm" not in tex
    assert tex.strip().startswith("$$") and tex.strip().endswith("$$")
    assert PAPER_SYMBOLS["tau"] == PAPER_SYMBOLS["psi"] == r"\psi"


def test_box_quadrature_exact_for_polynomial():
    f = lambda p: p[0] ** 2 * p[1] * p[2] ** 3
    val = gauss_legendre_box(f, np.zeros(3), np.ones(3), 4)
    assert abs(val - 1.0 / 24.0) < 1e-14


def test_half_ball_kernel_volume():
    r = 0.11
    one = lambda p: np.ones(p.shape[1:])
    _, vol = half_ball_kernel_integral(
        one, np.array([0.0, 0.5, 0.5]), r, np.zeros(3), np.ones(3), 40)
    assert abs(vol - 0.5 * math.pi * r**3 / 3.0) < 1e-8


def test_half_ball_kernel_interior_matches_full_ball():
    r = 0.1
    fx = lambda p: p[0]
    c = np.array([0.5, 0.5, 0.5])
    s, vol = half_ball_kernel_integral(fx, c, r, np.zeros(3), np.ones(3), 40)
    assert abs(vol - math.pi * r**3 / 3.0) < 1e-8
    assert abs(s - 0.5 * vol) < 1e-8
