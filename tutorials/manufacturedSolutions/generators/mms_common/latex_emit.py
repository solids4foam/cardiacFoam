"""LaTeX fragments in the paper's notation."""

import math

import sympy as sp

PAPER_SYMBOLS = {
    "Vm": r"V_m",
    "G": r"\boldsymbol\sigma",
    "M": r"\mathbf M",
    "tau": r"\psi",
    "psi": r"\psi",
    "cm": r"C_m",
    "chi": r"\chi",
    "c0": r"c_0",
    "beta": r"\beta",
    "u1": r"u_1",
    "u2": r"u_2",
    "u3": r"u_3",
}


def paper_symbol_names(symbols):
    """Map each symbol whose name is a PAPER_SYMBOLS key to its LaTeX."""
    return {s: PAPER_SYMBOLS[s.name] for s in symbols if s.name in PAPER_SYMBOLS}


def emit_latex_fragment(lhs, expr, symbol_names):
    """One displayed equation `lhs = expr` with symbols named by symbol_names."""
    body = sp.latex(expr, symbol_names=dict(symbol_names))
    return "$$\n%s = %s\n$$\n" % (lhs, body)


def format_value(value):
    """LaTeX of a parameter value: a string is kept, a sympy expression is
    printed, a number is printed in its shortest form (pi as \\pi)."""
    if isinstance(value, str):
        return value
    if isinstance(value, sp.MatrixBase) and 1 in value.shape:
        return r"\left(%s\right)" % ", ".join(sp.latex(v) for v in value)
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        return sp.latex(value)
    if isinstance(value, float) and abs(value - math.pi) < 1e-15:
        return r"\pi"
    text = "%.12g" % value
    if "e" in text:
        mantissa, exponent = text.split("e")
        text = r"%s\times10^{%d}" % (mantissa, int(exponent))
    return text


def emit_parameter_fragment(pairs, per_row=4):
    """One display of the `name = value` pairs, per_row to a row."""
    items = ["%s = %s" % (lhs, format_value(v)) for lhs, v in pairs]
    rows = [r",\quad ".join(items[i:i + per_row]) for i in range(0, len(items), per_row)]
    if len(rows) == 1:
        body = rows[0]
    else:
        body = "\\begin{gathered}\n" + ",\\\\\n".join(rows) + "\n\\end{gathered}"
    return "$$\n%s\n$$\n" % body
