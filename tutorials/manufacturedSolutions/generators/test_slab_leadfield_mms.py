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


# ---------------------------------------------------- C++ header regions ----

@pytest.mark.parametrize("name", sorted(slm.HEADER_PATHS))
def test_header_up_to_date(name):
    assert slm.header_matches(slm.HEADER_PATHS[name])


@pytest.mark.parametrize("name", sorted(slm.HEADER_PATHS))
def test_header_check_cli(name):
    assert slm.main(["--check", slm.HEADER_PATHS[name]]) == 0


def test_header_write_touches_only_regions(tmp_path):
    import mms_codegen as codegen
    for name, rel in slm.HEADER_PATHS.items():
        text = open(codegen.resolve(slm.__file__, rel)).read()
        stale = codegen.splice_regions(text, {n: ["// stale"] for n in slm.HEADERS[name](text)})
        assert stale != text
        dst = tmp_path/name
        dst.write_text(stale)
        assert slm.main(["--check", str(dst)]) == 1
        assert slm.main(["--write", str(dst)]) == 0
        assert dst.read_text() == text


def test_header_without_regions_is_refused(tmp_path):
    for name in slm.HEADER_PATHS:
        dst = tmp_path/name
        dst.write_text("namespace Foam {}\n")
        assert slm.main(["--check", str(dst)]) == 1
        assert slm.main(["--write", str(dst)]) == 1
        assert dst.read_text() == "namespace Foam {}\n"


def test_docs_reference_uses_the_generated_module():
    """docs/leadField/leadfield_reference.py takes the slab and layered
    potentials from reference/slab_leadfield.py, not from a transcription."""
    lf = _load(DOCS_REFERENCE, "leadfield_reference")
    ref = _reference()
    assert os.path.samefile(lf.SLAB.__file__, os.path.join(HERE, "reference", "slab_leadfield.py"))
    xs = np.concatenate([np.linspace(-1.0, 2.0, 31), [-0.5, 0.0, 1.0, 1.5]])
    largs = [lf.LAYERS[n] for n in ref.PARAMS["LAYER_NAMES"]]
    for tt in (0.0, 0.05):
        assert np.array_equal(lf.slab_phi_raw(xs, tt), ref.phi_slab(xs, tt, lf.K))
        assert np.array_equal(lf.layered_phi_raw(xs, tt),
                              ref.phi_layered(xs, tt, lf.K, lf.ALPHA, *largs))


def _parity_points():
    from parity.test_parity import sample_points
    p = sample_points()
    p[:, 0] = 3.0*p[:, 0] - 1.0
    edges = np.array([[v, 0.31, 0.62, 0.04] for v in (-1.0, -0.5, 0.0, 1.0, 1.5, 2.0)])
    return np.vstack([p, edges])


def _parity_or_skip():
    try:
        from parity.test_parity import compare, parity_available
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = parity_available()
    if not ok:
        pytest.skip(why)
    return compare


def test_parity_slab():
    compare = _parity_or_skip()
    ref = _reference()
    k = ref.PARAMS["K"]
    compare("leadFieldSlab", 1,
            lambda p: np.column_stack([ref.Vm(p[:, 0], p[:, 3]), ref.phi_slab(p[:, 0], p[:, 3], k)]),
            _parity_points(), rtol=1e-12, k="%.17g" % k)


def test_parity_layered():
    compare = _parity_or_skip()
    ref = _reference()
    k, a = ref.PARAMS["K"], ref.PARAMS["ALPHA"]
    largs = [ref.PARAMS["LAYERS"][n] for n in ref.PARAMS["LAYER_NAMES"]]
    keys = {n: "%.17g" % ref.PARAMS["LAYERS"][n] for n in ref.PARAMS["LAYER_NAMES"]}
    compare("leadFieldLayered", 1,
            lambda p: np.column_stack([ref.resistance(p[:, 0], *largs),
                                       ref.phi_layered(p[:, 0], p[:, 3], k, a, *largs)]),
            _parity_points(), rtol=1e-12, k="%.17g" % k, alpha="%.17g" % a, **keys)


def _case_value(path, pattern):
    with open(path) as fh:
        m = re.search(pattern, fh.read())
    return m.group(1)


@pytest.mark.skipif(not os.path.isdir(CASE), reason="tutorial case absent (cardiacFoam-only)")
def test_tutorial_matches():
    props = os.path.join(CASE, "constant", "electroProperties")
    # The set anisotropic: the xx entries of the ECG domain (monodomain heart)
    # and of bidomainSolverCoeffs
    tensor = r"\s*\[[^\]]*\]\s*\(\s*(\S+)"
    ecg = r"\becgDomains\s*\{[^}]*?"
    bidomain = r"\bbidomainSolverCoeffs\s*\{[^}]*?"
    gi = _case_value(props, ecg + r"conductivityIntracellular" + tensor)
    ge = _case_value(props, ecg + r"conductivityExtracellular" + tensor)
    assert (sp.Rational(gi), sp.Rational(ge)) == slm.MONODOMAIN_PAIR
    bi = _case_value(props, bidomain + r"conductivityIntracellular" + tensor)
    be = _case_value(props, bidomain + r"conductivityExtracellular" + tensor)
    assert (sp.Rational(bi), sp.Rational(be)) == slm.BIDOMAIN_PAIR
    pin = _case_value(props, r"\btorso\s*\{[^{}]*?phiERefPoint\s*\(([^)]*)\)")
    assert tuple(sp.Rational(v) for v in pin.split()) == slm.PHI_REF_POINT
    slab = os.path.join(CASE, "system", "setTorsoOrganConductivityFieldDict")
    assert sp.Rational(_case_value(slab, r"\bbath\s+(\S+);")) == slm.SIGMA_SLAB
    layered = slab + ".layered"
    assert sp.Rational(_case_value(layered, r"torsoInner\s+(\S+);")) == slm.SIGMA_INNER
    assert sp.Rational(_case_value(layered, r"torsoOuter\s+(\S+);")) == slm.SIGMA_OUTER
    allrun_path = os.path.join(CASE, "Allrun")
    allrun = open(allrun_path).read()
    assert sp.Rational(_case_value(allrun_path, r"\nlayerLeft=(\S+)")) == slm.X_LEFT
    assert sp.Rational(_case_value(allrun_path, r"\nlayerRight=(\S+)")) == slm.X_RIGHT
    assert "alpha=${alpha:-0.01}" in allrun
