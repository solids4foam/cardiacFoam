"""LaTeX fragments in the paper's notation."""

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
