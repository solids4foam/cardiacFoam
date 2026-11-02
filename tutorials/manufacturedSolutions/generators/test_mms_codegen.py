"""Marked-region splicing of hand-written C++ headers (mms_codegen)."""
import pytest
import sympy as sp

import mms_codegen as codegen

TEXT = """\
inline scalar f(const scalar x)
{
    if (x > 0)
    {
        // BEGIN GENERATED f_pos
        return old;
        // END GENERATED f_pos
    }
    // BEGIN GENERATED f_neg
    // END GENERATED f_neg
    return 0;   // hand-written
}
"""


def test_splice_replaces_only_region_bodies():
    out = codegen.splice_regions(TEXT, {"f_pos": ["return x;"], "f_neg": ["x = -x;", "y = 1;"]})
    assert out == TEXT.replace("        return old;\n", "        return x;\n").replace(
        "    // BEGIN GENERATED f_neg\n",
        "    // BEGIN GENERATED f_neg\n    x = -x;\n    y = 1;\n")


def test_splice_is_idempotent_and_leaves_other_regions():
    once = codegen.splice_regions(TEXT, {"f_pos": ["return x;"]})
    assert "    // BEGIN GENERATED f_neg\n    // END GENERATED f_neg\n" in once
    assert codegen.splice_regions(once, {"f_pos": ["return x;"]}) == once


def test_blank_generated_line_has_no_trailing_space():
    out = codegen.splice_regions(TEXT, {"f_neg": ["a;", "", "b;"]})
    assert "    a;\n\n    b;\n" in out


@pytest.mark.parametrize("text, message", [
    (TEXT, "missing"),
    (TEXT.replace("// END GENERATED f_pos", ""), "END"),
    (TEXT + "// BEGIN GENERATED f_pos\n// END GENERATED f_pos\n", "more than once"),
])
def test_splice_rejects_bad_markers(text, message):
    regions = {"f_pos": ["return x;"]}
    if message == "missing":
        regions = {"f_other": ["return x;"]}
    with pytest.raises(ValueError, match=message):
        codegen.splice_regions(text, regions)


def test_long_assignment_wraps_at_top_level_terms():
    a, b, c = sp.symbols("alpha_long_name beta_long_name gamma_long_name", real=True)
    lines = codegen.assignment("value", sp.sin(a)*sp.cos(b) + a*b*c + sp.exp(c)*a + b**3*c, 4)
    assert len(lines) > 2 and all(4 + len(line) <= 79 for line in lines)
    joined = " ".join(l.strip() for l in lines)
    assert joined.startswith("value =") and joined.endswith(";")


def test_print_pi_as_local():
    x = sp.Symbol("x", real=True)
    assert codegen.cxx_local(sp.cos(sp.pi*x)) == "cos(pi*x)"


REGION_GENERATORS = ("fda_monodomain_mms", "pseudo_ecg_mms", "rotated_monodomain_mms",
                     "cartesian_eikonal_mms", "eikonal_ecg_mms")
REGION_HEADERS = ("manufacturedFDAReference.H",
                  "manufacturedAnisotropicMonodomainReference.H",
                  "manufacturedEikonalReference.H")


@pytest.mark.parametrize("header", REGION_HEADERS)
def test_every_marker_owned_by_exactly_one_generator(header):
    import importlib
    owners, path = {}, None
    for name in REGION_GENERATORS:
        mod = importlib.import_module(name)
        if header in mod.HEADERS:
            path = codegen.resolve(mod.__file__, mod.HEADER_PATHS[header])
            text = open(path).read()
            for region in mod.HEADERS[header](text):
                owners.setdefault(region, []).append(name)
    assert path is not None
    twice = {r: o for r, o in owners.items() if len(o) > 1}
    assert not twice, "regions owned twice: %r" % twice
    marked = set(codegen._regions_in(open(path).read().split("\n")))
    assert marked == set(owners), (
        "marked but unowned: %r; owned but unmarked: %r"
        % (sorted(marked - set(owners)), sorted(set(owners) - marked)))


# Every generator in the release: --emit and --self-check, so the README can
# document one command per case; --check/--write only where a header is owned.
CLI_GENERATORS = ("fda_monodomain_mms", "rotated_monodomain_mms", "pseudo_ecg_mms",
                  "graph_monodomain_mms", "coupled_1d3d_mms", "cartesian_eikonal_mms",
                  "eikonal_ecg_mms", "shell_monodomain_mms", "shell_eikonal_mms",
                  "template_mms")


def _stubbed(monkeypatch, name, fail=False):
    import importlib
    mod = importlib.import_module(name)
    calls = []

    def checks(*args, **kwargs):
        calls.append("checks")
        if fail:
            raise AssertionError("stub failure")

    def outputs(*args, **kwargs):
        calls.append("outputs")
        return []

    monkeypatch.setattr(mod, "run_self_checks", checks)
    monkeypatch.setattr(mod, "write_outputs", outputs)
    return mod, calls


@pytest.mark.parametrize("name", CLI_GENERATORS)
def test_cli_self_check_and_emit(monkeypatch, name):
    mod, calls = _stubbed(monkeypatch, name)
    assert mod.main(["--self-check"]) == 0 and calls == ["checks"]
    del calls[:]
    assert mod.main(["--emit"]) == 0 and calls == ["checks", "outputs"]


@pytest.mark.parametrize("name", CLI_GENERATORS)
def test_cli_self_check_failure_exits_nonzero(monkeypatch, name):
    mod, calls = _stubbed(monkeypatch, name, fail=True)
    assert mod.main(["--self-check"]) != 0
    assert mod.main(["--emit"]) != 0 and "outputs" not in calls


@pytest.mark.parametrize("name", ("graph_monodomain_mms", "coupled_1d3d_mms", "template_mms"))
def test_cli_without_header_has_no_check_path(monkeypatch, name):
    mod, _ = _stubbed(monkeypatch, name)
    with pytest.raises(SystemExit):
        mod.main(["--check", "x.H"])


def test_rationals_print_as_floating_point():
    x = sp.Symbol("x", real=True)
    assert codegen.cxx_local(sp.Rational(1, 2)*x) == "(1.0/2.0)*x"
    assert "/" not in codegen.cxx_local(x**3).replace("1.0/", "")
