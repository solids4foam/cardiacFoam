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


def test_rationals_print_as_floating_point():
    x = sp.Symbol("x", real=True)
    assert codegen.cxx_local(sp.Rational(1, 2)*x) == "(1.0/2.0)*x"
    assert "/" not in codegen.cxx_local(x**3).replace("1.0/", "")
