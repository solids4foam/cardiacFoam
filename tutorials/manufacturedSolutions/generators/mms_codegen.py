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


# ---------------------------------------------------------------------------
# Generated regions inside hand-written headers.
#
# A hand-written header marks each generated formula with
#
#     // BEGIN GENERATED <name>
#     ...
#     // END GENERATED <name>
#
# A generator owns a set of names; it rewrites the lines between its markers,
# indented as the BEGIN marker, and leaves every other byte of the file as it
# is.  Region lines are given without that indentation.
# ---------------------------------------------------------------------------

LINE_WIDTH = 79

_BEGIN = re.compile(r"^([ \t]*)// BEGIN GENERATED (\S+)[ \t]*$")
_END = re.compile(r"^[ \t]*// END GENERATED (\S+)[ \t]*$")


class _LocalPiPrinter(_CxxPrinter):
    """As _CxxPrinter, with pi printed as a local `pi` and the square of a
    compound base as sqr(base)."""

    def _print_Pi(self, expr):
        return "pi"

    def _print_Pow(self, expr):
        if expr.exp == 2 and not expr.base.is_Symbol:
            return "sqr(%s)" % self._print(expr.base)
        return super()._print_Pow(expr)


_LOCAL_PI_PRINTER = _LocalPiPrinter()


def cxx_local(expr):
    """C++ of expr with pi printed as `pi` (the function declares
    `const scalar pi = constant::mathematical::pi;`).  Raises on an integer
    division, which C++ would truncate."""
    code = _LOCAL_PI_PRINTER.doprint(expr)
    if re.search(r"(?<![\w.])\d+/\d+(?![\w.])", code):
        raise RuntimeError("integer division in emitted code: " + code)
    return code


def _top_level_split(code, separators):
    """Split code at the separators that sit outside parentheses; returns
    [(separator or '', term), ...]."""
    out, depth, start, i, sep = [], 0, 0, 0, ""
    while i < len(code):
        ch = code[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif depth == 0:
            for s in separators:
                if code.startswith(s, i) and i > 0:
                    out.append((sep, code[start:i]))
                    sep, start = s.strip(), i + len(s)
                    i += len(s) - 1
                    break
        i += 1
    out.append((sep, code[start:]))
    return out


def _block(code, col, width):
    """[(column, text)] of code starting at column col (>= 4): broken at
    top-level + and - (operator two columns left of the terms), else at
    top-level * and / (operator one column left), until each line fits width.  A
    factor that alone is too long stays on one line.  Only whitespace is
    inserted, so the C++ is unchanged."""
    if col + len(code) <= width:
        return [(col, code)]
    for separators, back in (((" + ", " - "), 2), (("*", "/"), 1)):
        parts = _top_level_split(code, separators)
        if len(parts) > 1:
            break
    else:
        return [(col, code)]
    out = []
    for k, (sep, term) in enumerate(parts):
        sub = _block(term, col, width)
        if k:
            head = sep + " " if back == 2 else sep
            sub[0] = (col - back, head + sub[0][1])
        out += sub
    return out


def wrap_expression(code, indent):
    """Lines (indented relative to a statement at column indent) of code
    placed under that statement."""
    return [" "*c + text for c, text in _block(code, 4, LINE_WIDTH - indent)]


def assignment(lhs, expr, indent, decl=""):
    """Lines of `decl lhs = expr;` (pi printed as local), wrapped for a
    region at column indent."""
    code = expr if isinstance(expr, str) else cxx_local(expr)
    head = ("%s %s" % (decl, lhs)).strip()
    one = "%s = %s;" % (head, code)
    if indent + len(one) <= LINE_WIDTH:
        return [one]
    body = wrap_expression(code, indent)
    body[-1] += ";"
    return ["%s =" % head] + body


def statement(keyword, expr, indent):
    """Lines of `keyword expr;` (e.g. return), wrapped like assignment."""
    code = expr if isinstance(expr, str) else cxx_local(expr)
    one = "%s %s;" % (keyword, code)
    if indent + len(one) <= LINE_WIDTH:
        return [one]
    body = wrap_expression(code, indent)
    body[-1] += ";"
    return [keyword] + body


def call(name, args, indent, keyword="return"):
    """Lines of `keyword name(arg, ...);`, one argument per line when long."""
    codes = [a if isinstance(a, str) else cxx_local(a) for a in args]
    one = "%s %s(%s);" % (keyword, name, ", ".join(codes))
    if indent + len(one) <= LINE_WIDTH:
        return [one]
    lines = ["%s %s" % (keyword, name), "("]
    for i, c in enumerate(codes):
        sub = wrap_expression(c, indent)
        sub[-1] += "," if i < len(codes) - 1 else ""
        lines += sub
    lines.append(");")
    return lines


def cse_lines(exprs, indent, prefix="w"):
    """(declaration lines, reduced exprs) of the common subexpressions of
    exprs, named prefix0, prefix1, ..."""
    repl, red = sp.cse([sp.sympify(e) for e in exprs],
                       symbols=sp.numbered_symbols(prefix, real=True), order="canonical")
    lines = []
    for sym, val in repl:
        lines += assignment(str(sym), val, indent, decl="const scalar")
    return lines, red


def _regions_in(lines):
    """{name: (begin line index, end line index, indentation)}."""
    found, open_name = {}, None
    for i, line in enumerate(lines):
        b, e = _BEGIN.match(line), _END.match(line)
        if b:
            if open_name is not None:
                raise ValueError("BEGIN GENERATED %s has no END GENERATED before line %d"
                                 % (open_name[0], i + 1))
            name = b.group(2)
            if name in found:
                raise ValueError("region %s is marked more than once" % name)
            open_name = (name, i, b.group(1))
        elif e:
            if open_name is None or e.group(1) != open_name[0]:
                raise ValueError("END GENERATED %s at line %d does not close the open region"
                                 % (e.group(1), i + 1))
            found[open_name[0]] = (open_name[1], i, open_name[2])
            open_name = None
    if open_name is not None:
        raise ValueError("BEGIN GENERATED %s has no END GENERATED" % open_name[0])
    return found


def splice_regions(text, regions):
    """text with the body of each region in regions ({name: [lines]})
    replaced; raises ValueError when a region is missing or the markers are
    malformed."""
    lines = text.split("\n")
    found = _regions_in(lines)
    missing = sorted(set(regions) - set(found))
    if missing:
        raise ValueError("regions missing from the header: %s" % ", ".join(missing))
    out, i = [], 0
    for name in sorted(regions, key=lambda n: found[n][0]):
        begin, end, pad = found[name]
        out += lines[i:begin + 1]
        out += [(pad + s) if s else "" for s in regions[name]]
        i = end
    out += lines[i:]
    return "\n".join(out)


def region_indent(text, name):
    """Indentation (columns) of the BEGIN marker of region name."""
    found = _regions_in(text.split("\n"))
    if name not in found:
        raise ValueError("regions missing from the header: %s" % name)
    return len(found[name][2])


def write_regions(module_file, path, regions_for):
    """Rewrite the regions of the header at path; regions_for(text) returns
    {name: [lines]}.  Returns the path written."""
    path = resolve(module_file, path)
    with open(path) as fh:
        text = fh.read()
    new = splice_regions(text, regions_for(text))
    if new != text:
        with open(path, "w") as fh:
            fh.write(new)
    return path


def regions_match(module_file, path, regions_for):
    """True when the regions of the header at path equal the generator's."""
    path = resolve(module_file, path)
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        text = fh.read()
    try:
        return splice_regions(text, regions_for(text)) == text
    except ValueError:
        return False


def generator_main(argv, doc, module_file, run_self_checks, write_outputs, headers):
    """Command line of a generator with emitted files and header regions:
    --emit, --self-check, --write PATH, --check PATH.  --emit and --write run
    the self-checks first."""
    import argparse
    ap = argparse.ArgumentParser(description=doc.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--emit", action="store_true",
                     help="write reference/*.py and latex/*.tex")
    grp.add_argument("--self-check", action="store_true", help="run the self-checks")
    grp.add_argument("--write", metavar="PATH",
                     help="rewrite the generated regions of a C++ header")
    grp.add_argument("--check", metavar="PATH",
                     help="compare the generated regions of a C++ header")
    args = ap.parse_args(argv)
    if args.check:
        return regions_cli(module_file, args.check, headers, write=False)
    try:
        run_self_checks()
    except AssertionError as err:
        print("self-check failed: %s" % err)
        return 1
    if args.write:
        return regions_cli(module_file, args.write, headers, write=True)
    if args.emit:
        for rel in write_outputs():
            print("wrote %s" % rel)
        return 0
    print("self-checks passed")
    return 0


def regions_cli(module_file, path, headers, write):
    """--write / --check of a hand-written header with generated regions.

    headers maps a header file name to its regions_for(text) function; the
    file name of path selects it.  Returns the exit status."""
    name = os.path.basename(path)
    if name not in headers:
        print("%s holds no regions of this generator; known headers: %s"
              % (name, ", ".join(sorted(headers))))
        return 2
    full = resolve(module_file, path)
    if write:
        try:
            write_regions(module_file, path, headers[name])
        except ValueError as err:
            print("%s: %s" % (full, err))
            return 1
        print("wrote the generated regions of %s" % full)
        return 0
    if regions_match(module_file, path, headers[name]):
        print("%s is up to date" % full)
        return 0
    print("%s differs from the generator output (or lacks its regions)" % full)
    return 1
