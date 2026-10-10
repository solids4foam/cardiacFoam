"""C++ emission shared by the transmission manufactured-solution generators.

A geometry module (transmission_mms, shell_mms) provides the symbolic fields
and these attributes:

    ARGS, P, PARAM_ORDER, SCALAR_PARAMS, TAIL_PARAMS, ARRAY_PARAMS,
    ANGLE_PARAMS, PREAMBLE, PREAMBLE_NEEDS, SOLUTIONS, HEART_ZONE, ZONE_NAMES,
    ZONES_OF_KIND, SIGNATURES, RETURN_TYPE, _expr(solution, zone, name),
    printed_check_context()  -> {argument name: float}

and the functions below turn them into the reference header.
"""

import math
import os
import re

import sympy as sp
from sympy.printing.c import C99CodePrinter


class _CxxPrinter(C99CodePrinter):
    """C++ for OpenFOAM scalars: no integer-exponent pow(), pi from Foam."""

    def _print_Pi(self, expr):
        return "constant::mathematical::pi"

    def _print_Pow(self, expr):
        # Every result is a single parenthesised or call operand, so it is
        # safe after '/' in the Mul printer.
        b, e = expr.base, expr.exp
        bs = self._print(b) if b.is_Symbol else "(%s)" % self._print(b)
        if e.is_Integer and 1 <= abs(int(e)) <= 4:
            prod = "*".join([bs]*abs(int(e)))
            return "(%s)" % prod if e > 0 else "(1.0/(%s))" % prod
        if e == sp.Rational(1, 2):
            return "sqrt(%s)" % self._print(b)
        if e == -sp.Rational(1, 2):
            return "(1.0/sqrt(%s))" % self._print(b)
        return "pow(%s, %s)" % (self._print(b), self._print(sp.Float(e, 17)))


_PRINTER = _CxxPrinter()


def cxx(expr):
    return _PRINTER.doprint(expr)


def as_list(e):
    return list(e) if isinstance(e, (list, tuple, sp.MatrixBase)) else [e]


def g17(v):
    return "%.17g" % v


def check_printed(geo, exprs, repl, red):
    """Evaluate the printed C++ arithmetic as Python at one point and compare
    with sympy, so a printer precedence error cannot reach the header."""
    env = {"sin": math.sin, "cos": math.cos, "exp": math.exp, "sqrt": math.sqrt,
           "pow": math.pow}
    ctx = geo.printed_check_context()
    env.update(ctx)

    def py(code):
        if re.search(r"(?<![\w.])\d+/\d+(?![\w.])", code):
            raise RuntimeError("integer division in emitted code: " + code)
        return eval(code.replace("constant::mathematical::pi", "%r" % math.pi), {}, env)

    for sym, val in repl:
        env[str(sym)] = py(cxx(val))
    ref = sp.lambdify(geo.ARGS, exprs, modules="math")(*[ctx[str(a)] for a in geo.ARGS])
    for e, r in zip(red, ref):
        v = py(cxx(e))
        if abs(v - r) > 1e-12*max(1.0, abs(r)):
            raise RuntimeError("printed C++ differs from sympy: %r vs %r" % (v, r))


def branch_code(geo, exprs, vkind, indent):
    exprs = [sp.sympify(e) for e in exprs]
    repl, red = sp.cse(exprs, symbols=sp.numbered_symbols("w"), order="canonical")
    check_printed(geo, exprs, repl, red)
    pad = " "*indent
    lines = ["%sconst scalar %s = %s;" % (pad, cxx(sym), cxx(val)) for sym, val in repl]
    r = [cxx(e) for e in red]
    if vkind == "scalar":
        lines.append("%sreturn %s;" % (pad, r[0]))
    elif vkind == "vector":
        lines.append("%sreturn vector(%s, %s, %s);" % (pad, *r))
    elif vkind == "symmTensor":
        lines.append("%sreturn symmTensor(%s, %s, %s, %s, %s, %s);" % (pad, *r))
    else:
        lines.append("%sreturn n.x()*(%s) + n.y()*(%s) + n.z()*(%s);" % (pad, *r))
    return lines


def param_access(geo, name):
    for arr, _ in geo.ARRAY_PARAMS:
        if name.startswith(arr) and name[len(arr):].isdigit():
            return "p.%s[%s]" % (arr, name[len(arr):])
    return "p." + name


def emit_function(geo, idx, cname, kind, ename, vkind):
    zones = geo.ZONES_OF_KIND[kind]
    branches = {}
    used = set()
    for sol in geo.SOLUTIONS:
        for zn in zones:
            ex = as_list(geo._expr(sol, zn, ename))
            branches[(sol, zn)] = ex
            for e in ex:
                used |= sp.sympify(e).free_symbols
    for sym, needs in geo.PREAMBLE_NEEDS.items():
        if sym in used:
            used |= set(needs)
    pre = [(sym, deps, code) for sym, deps, code in geo.PREAMBLE if sym in used]
    params = [n for n in geo.PARAM_ORDER
              if geo.P[n] in used or any(n in deps for _, deps, _ in pre)]
    out = ["inline %s %s(%s)" % (geo.RETURN_TYPE[vkind], cname, geo.SIGNATURES[kind]), "{"]
    out += ["    const scalar %s = %s;" % (n, param_access(geo, n)) for n in params]
    out += ["    " + code for _, _, code in pre]
    if kind in ("zone", "bath"):
        out.append("    checkZone(zone, %s);" % ("true" if kind == "zone" else "false"))
    out.append("")

    def zone_block(sol, indent):
        lines = []
        if len(zones) == 1 or kind == "tensor":
            return branch_code(geo, branches[(sol, zones[0])], vkind, indent)
        pad = " "*indent
        for zn in zones:
            if zn == geo.HEART_ZONE:
                continue
            lines.append("%sif (zone == %s)" % (pad, geo.ZONE_NAMES[zn]))
            lines.append(pad + "{")
            lines += branch_code(geo, branches[(sol, zn)], vkind, indent + 4)
            lines.append(pad + "}")
        if geo.HEART_ZONE in zones:
            lines += branch_code(geo, branches[(sol, geo.HEART_ZONE)], vkind, indent)
        else:
            lines.append("%sreturn 0;" % pad)
        return lines

    for i, sol in enumerate(geo.SOLUTIONS[:-1]):
        out.append("    if (p.solution == %d)" % i)
        out.append("    {")
        out += zone_block(sol, 8)
        out.append("    }")
        out.append("")
    out += zone_block(geo.SOLUTIONS[-1], 4)
    out.append("}")
    return out


def params_function(geo, name, params, solution):
    v = geo.param_values(params)
    out = ["inline Params %s()" % name, "{", "    Params p;"]
    for n in geo.SCALAR_PARAMS + geo.TAIL_PARAMS:
        if n in geo.ANGLE_PARAMS:
            out.append("    p.%s = degToRad(%r);" % (n, float(params[n])))
        else:
            out.append("    p.%s = %r;" % (n, v[n]))
    for arr, cnt in geo.ARRAY_PARAMS:
        for i in range(cnt):
            out.append("    p.%s[%d] = %r;" % (arr, i, float(params[arr][i])))
    out.append("    p.solution = %d;" % solution)
    out += ["    return p;", "}"]
    return out


INCLUDES = ['#include "error.H"', '#include "label.H"', '#include "mathematicalConstants.H"',
            '#include "point.H"', '#include "scalar.H"', '#include "symmTensor.H"',
            '#include "unitConversion.H"', '#include "vector.H"', ""]


def check_zone_function(first, last, heart):
    return ["inline void checkZone(const label zone, const bool heartAllowed)", "{",
            "    if (zone < %s || zone > %s || (!heartAllowed && zone == %s))"
            % (first, last, heart),
            "    {", "        FatalErrorInFunction",
            "            << \"Zone \" << zone << \" is not valid here\" << exit(FatalError);",
            "    }", "}", ""]


def self_test_table(rows):
    L = ["struct SelfTestSample", "{", "    label fn;", "    label param;",
         "    label solution;", "    label zone;", "    label component;",
         "    scalar x, y, z, t, value;", "};", ""]
    L.append("inline const SelfTestSample selfTestSamples[] =")
    L.append("{")
    for r in rows:
        L.append("    {%d, %d, %d, %d, %d, %s, %s, %s, %s, %s}," % (
            r[0], r[1], r[2], r[3], r[4], g17(r[5]), g17(r[6]), g17(r[7]), g17(r[8]),
            g17(r[9])))
    L.append("};")
    L += ["", "inline const label nSelfTestSamples = %d;" % len(rows), ""]
    return L


def resolve(module_file, path):
    """Relative paths are taken against the generator directory, for --write
    and --check alike, so the result does not depend on the working directory."""
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(module_file)), path))


def header_matches(module_file, generate, path, text=None):
    path = resolve(module_file, path)
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        current = fh.read()
    return current == (text if text is not None else generate())
