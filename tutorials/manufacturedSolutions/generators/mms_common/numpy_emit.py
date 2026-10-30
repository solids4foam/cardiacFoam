"""Emission of NumPy modules from sympy expressions.

The emitted module imports only numpy and its functions broadcast over arrays.
"""

import sympy as sp
from sympy.printing.numpy import NumPyPrinter


def _printer():
    return NumPyPrinter({"inline": True})


def _function_source(name, args, expr):
    """Source of one vectorised function returning a scalar-field array or a
    stacked array of shape (rows, cols..., points) for a Matrix (column vectors drop the column axis)."""
    printer = _printer()
    arg_names = [printer.doprint(a) for a in args]
    is_matrix = isinstance(expr, sp.MatrixBase)
    entries = list(expr) if is_matrix else [expr]
    temps, reduced = sp.cse(entries, symbols=sp.numbered_symbols("_t"))
    lines = ["def %s(%s):" % (name, ", ".join(arg_names))]
    if arg_names:
        lines.append("    _ref = numpy.zeros(numpy.broadcast(%s).shape)" % ", ".join(arg_names))
    else:
        lines.append("    _ref = 0.0")
    for sym, sub in temps:
        lines.append("    %s = %s" % (sym.name, printer.doprint(sub)))
    outs = ["(%s) + _ref" % printer.doprint(e) for e in reduced]
    if is_matrix:
        lines.append("    _out = numpy.stack([%s])" % ", ".join(outs))
        rows, cols = expr.shape
        lead = (rows,) if cols == 1 else (rows, cols)
        lines.append("    return _out.reshape(%r + numpy.shape(_ref))" % (lead,))
    else:
        lines.append("    return %s" % outs[0])
    return "\n".join(lines)


def emit_numpy_module(module_name, functions, params, doc):
    """Source text of a numpy-only module.

    functions maps a function name to (arguments, expression); a Matrix
    expression returns an array of shape (rows, [cols,] ...).  params become
    the PARAMS dict.
    """
    parts = ['"""%s"""' % doc.replace('"""', "'''"), "",
             "import numpy", "", "# module: %s" % module_name, "",
             "PARAMS = %r" % (dict(params),), ""]
    for name, (args, expr) in functions.items():
        parts.extend(["", _function_source(name, args, expr), ""])
    return "\n".join(parts).rstrip("\n") + "\n"
