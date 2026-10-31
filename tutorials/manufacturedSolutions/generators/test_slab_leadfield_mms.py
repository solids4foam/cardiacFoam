import importlib.util
import os
import re
import sys

import numpy as np
import pytest
import sympy as sp

import slab_leadfield_mms as slm

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "..", "leadFieldECG")
DOCS_REFERENCE = os.path.join(HERE, "docs", "leadField", "leadfield_reference.py")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _reference():
    return _load(os.path.join(HERE, "reference", "slab_leadfield.py"), "slab_leadfield_ref")


def test_self_checks_and_cli():
    assert slm.main(["--self-check"]) == 0


def test_heart_equation_and_continuity():
    gi, ge = sp.symbols("gi ge", positive=True)
    heart = slm.slab_branches()[1][0].subs(slm.k, gi/(gi + ge))
    assert sp.simplify((gi + ge)*sp.diff(heart, slm.x, 2) + gi*sp.diff(slm.vm(), slm.x, 2)) == 0
    f = slm.fields()
    layers = slm.case_layers()
    R = sp.lambdify(slm.x, f["resistance"].subs(layers), "numpy")
    for xi in (-0.5, 0.0, 1.0, 1.5):
        assert abs(R(xi - 1e-12) - R(xi + 1e-12)) < 1e-9


def test_emitted_files_up_to_date():
    assert slm.outputs_match()


def test_emitted_files_regenerate(tmp_path):
    slm.write_outputs(str(tmp_path))
    for rel, text in slm.render_outputs().items():
        assert (tmp_path/rel).read_text() == text


def test_emitted_numpy_matches_sympy():
    ref = _reference()
    f = slm.fields()
    layers = slm.case_layers()
    xs = np.linspace(-1.0, 2.0, 61) + 1e-3
    for tt in (0.0, 0.05):
        k = ref.PARAMS["K"]
        sub = {slm.t: tt, slm.k: k, slm.alpha: ref.PARAMS["ALPHA"]}
        sub.update(layers)
        want = [float(f["phiLayered"].subs(sub).subs(slm.x, v)) for v in xs]
        largs = [ref.PARAMS["LAYERS"][n] for n in ref.PARAMS["LAYER_NAMES"]]
        got = ref.phi_layered(xs, tt, k, ref.PARAMS["ALPHA"], *largs)
        np.testing.assert_allclose(got, want, rtol=1e-13, atol=1e-15)


def test_numpy_matches_the_cxx_transcription():
    """reference/slab_leadfield.py against docs/leadField/leadfield_reference.py,
    which follows manufacturedLeadFieldECGReference.H branch by branch."""
    ref = _reference()
    lf = _load(DOCS_REFERENCE, "leadfield_reference")
    assert abs(ref.PARAMS["K"] - lf.K) < 1e-15
    assert ref.PARAMS["ALPHA"] == lf.ALPHA
    largs = [ref.PARAMS["LAYERS"][n] for n in ref.PARAMS["LAYER_NAMES"]]
    for n in ref.PARAMS["LAYER_NAMES"]:
        assert abs(ref.PARAMS["LAYERS"][n] - lf.LAYERS[n]) <= 1e-15*abs(lf.LAYERS[n]), n
    xs = np.concatenate([np.linspace(-1.0, 2.0, 301), [-0.5, 0.0, 1.0, 1.5]])
    for tt in (0.0, 0.013, 0.05):
        np.testing.assert_allclose(ref.Vm(xs, tt), lf.slab_vm(xs, tt), rtol=0, atol=1e-15)
        np.testing.assert_allclose(ref.phi_slab(xs, tt, lf.K), lf.slab_phi_raw(xs, tt),
                                   rtol=0, atol=1e-15)
        np.testing.assert_allclose(ref.resistance(xs, *largs), lf.layered_resistance(xs),
                                   rtol=1e-15, atol=1e-14)
        np.testing.assert_allclose(ref.phi_layered(xs, tt, lf.K, lf.ALPHA, *largs),
                                   lf.layered_phi_raw(xs, tt), rtol=1e-15, atol=1e-15)


def _case_value(path, pattern):
    with open(path) as fh:
        m = re.search(pattern, fh.read())
    return m.group(1)


@pytest.mark.skipif(not os.path.isdir(CASE), reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_matches():
    props = os.path.join(CASE, "constant", "electroProperties")
    gi = _case_value(props, r"_monodomainGi\s*\(\s*(\S+)")
    ge = _case_value(props, r"_monodomainGe\s*\(\s*(\S+)")
    assert (sp.Rational(gi), sp.Rational(ge)) == slm.MONODOMAIN_PAIR
    bi = _case_value(props, r"_bidomainGi\s*\(\s*(\S+)")
    be = _case_value(props, r"_bidomainGe\s*\(\s*(\S+)")
    assert (sp.Rational(bi), sp.Rational(be)) == slm.BIDOMAIN_PAIR
    pin = _case_value(props, r"_leadFieldTorso\s*\{[^{}]*?phiERefPoint\s*\(([^)]*)\)")
    assert tuple(sp.Rational(v) for v in pin.split()) == slm.PHI_REF_POINT
    slab = os.path.join(CASE, "system", "setTorsoOrganConductivityFieldDict")
    assert sp.Rational(_case_value(slab, r"\bbath\s+(\S+);")) == slm.SIGMA_SLAB
    layered = slab + ".layered"
    assert sp.Rational(_case_value(layered, r"torsoInner\s+(\S+);")) == slm.SIGMA_INNER
    assert sp.Rational(_case_value(layered, r"torsoOuter\s+(\S+);")) == slm.SIGMA_OUTER
    allrun = open(os.path.join(CASE, "Allrun")).read()
    assert "layerInterfaces ( -0.5 1.5 );" in allrun
    assert "alpha=${alpha:-0.01}" in allrun
