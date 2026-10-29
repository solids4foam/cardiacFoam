import importlib.util
import os
import sys

import numpy as np
import pytest
import sympy as sp

import template_mms as tpl

x, y, t = tpl.x, tpl.y, tpl.t


def test_residual_vanishes():
    assert tpl.residual() == 0


def test_source_is_chi_cm_dVdt_minus_div_sigma_grad_V():
    f = tpl.fields()
    V = f["Vm"]
    div = tpl.SIGMA[0, 0]*sp.diff(V, x, 2) + tpl.SIGMA[1, 1]*sp.diff(V, y, 2)
    assert sp.simplify(f["S"] - (tpl.CHI*tpl.CM*sp.diff(V, t) - div)) == 0


def test_conormal_flux_vanishes_on_every_side():
    fluxes = tpl.conormal_fluxes()
    assert sorted(fluxes) == [("x", 0), ("x", 1), ("y", 0), ("y", 1)]
    assert all(v == 0 for v in fluxes.values())


def test_source_matches_finite_differences():
    assert tpl.max_source_fd_error() < tpl.FD_RTOL


def test_printed_source_equals_source():
    S = tpl.printed_source(tpl.CHI, tpl.CM, tpl.SIGMA[0, 0], tpl.SIGMA[1, 1])
    assert sp.simplify(S - tpl.fields()["S"]) == 0


def test_self_checks_catch_a_source_without_diffusion(monkeypatch):
    V = tpl.exact_field()
    monkeypatch.setattr(tpl, "fields",
                        lambda *a, **k: {"Vm": V, "S": tpl.CHI*tpl.CM*sp.diff(V, t)})
    with pytest.raises(AssertionError):
        tpl.run_self_checks()


def test_self_checks_and_cli():
    assert tpl.main(["--self-check"]) == 0


def test_cli_has_no_header_options():
    with pytest.raises(SystemExit):
        tpl.main(["--check", "x.H"])


def test_emitted_files_up_to_date():
    assert tpl.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    tpl.write_outputs(str(tmp_path))
    for rel, text in tpl.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def _reference(monkeypatch):
    monkeypatch.setitem(sys.modules, "sympy", None)
    path = os.path.join(tpl.HERE, "reference", "template.py")
    spec = importlib.util.spec_from_file_location("template_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_emitted_module_is_numpy_only(monkeypatch):
    with open(os.path.join(tpl.HERE, "reference", "template.py")) as fh:
        assert "sympy" not in fh.read().split('"""')[-1]
    ref = _reference(monkeypatch)
    assert ref.PARAMS == {"CHI": 3.0, "CM": 2.0, "SIGMA": [0.13, 0.0176]}


def test_emitted_numpy_matches_sympy(monkeypatch):
    f = tpl.fields()
    fns = {n: sp.lambdify((x, y, t), e, "numpy") for n, e in f.items()}
    ref = _reference(monkeypatch)
    pts = np.random.default_rng(7).uniform(0, 1, (6, 3))
    for name, fn in fns.items():
        np.testing.assert_allclose(getattr(ref, name)(*pts.T), fn(*pts.T),
                                   rtol=1e-13, atol=1e-15)
