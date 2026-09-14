#!/usr/bin/env python3
"""Smoke assertions for the static immersed-cylinder regression case.

These are health checks on the immersed-boundary machinery, not a validation
of the flow. They are meant to catch a solver that diverges, produces NaN,
silently stops forcing the immersed body, or loses mass. They deliberately do
not check drag: this case resolves the cylinder with only a handful of cells
(see README.md), so no drag coefficient from it would be meaningful.

Limits are generous bands around values measured on OpenFOAM v2412 and v2512
(macOS/clang) and are documented alongside each check.

Standard library only.
"""

import math
import re
import sys

END_TIME = "0.02"
LOG = "log.pimpleHFDIBFoam"

failures = []
checks = 0


def check(name, ok, detail):
    global checks
    checks += 1
    if ok:
        print(f"  ok    {name}: {detail}")
    else:
        print(f"  FAIL  {name}: {detail}")
        failures.append(name)


def read_internal(path):
    """Return (kind, values) for the internalField of an ascii OpenFOAM field."""
    with open(path) as fh:
        txt = fh.read()
    m = re.search(
        r"internalField\s+nonuniform\s+List<(\w+)>\s*\n(\d+)\s*\n\(\n(.*?)\n\)\s*;",
        txt,
        re.S,
    )
    if m:
        kind, body = m.group(1), m.group(3)
        if kind == "scalar":
            return kind, [float(x) for x in body.split()]
        return kind, [
            tuple(float(v) for v in t.split()) for t in re.findall(r"\(([^)]*)\)", body)
        ]
    m = re.search(r"internalField\s+uniform\s+(.+?);", txt, re.S)
    if m:
        raw = m.group(1).strip()
        if raw.startswith("("):
            return "vector", [tuple(float(v) for v in raw.strip("()").split())]
        return "scalar", [float(raw)]
    raise SystemExit(f"could not parse internalField from {path}")


def mags(vectors):
    return [math.sqrt(sum(c * c for c in v)) for v in vectors]


try:
    with open(LOG) as fh:
        log = fh.read()
except OSError as exc:
    print(f"FAIL: cannot read {LOG}: {exc}")
    sys.exit(1)

# 1. The run reached its end time and terminated normally.
check(
    "reached end time",
    re.search(rf"^Time = {re.escape(END_TIME)}\b", log, re.M) is not None,
    f"'Time = {END_TIME}' present in {LOG}",
)
check("clean termination", log.rstrip().endswith("End"), "log ends with 'End'")

# 2. Exactly one immersed body was active, at every report.
counts = [int(n) for n in re.findall(r"Active IB listZize\s*:\s*(\d+)", log)]
check(
    "immersed body count",
    len(counts) > 0 and set(counts) == {1},
    f"{len(counts)} report(s), values {sorted(set(counts)) or 'none'}, expected all 1",
)

# 3. No solver-reported failure. Scanned line by line so that OpenFOAM's
# startup banner, "trapFpe: Floating point exception trapping enabled", is not
# mistaken for an actual exception. Field values are checked numerically below;
# this only catches what the solver itself reported.
bad_lines = [
    line
    for line in log.splitlines()
    if "trapFpe" not in line
    and re.search(r"\b(nan|inf)\b|Floating point exception", line, re.I)
]
check(
    "no floating point error in log",
    not bad_lines,
    "no nan/inf/FPE reported" if not bad_lines else f"first offending line: {bad_lines[0][:90]!r}",
)

# 4. Mass conservation stayed bounded.
cum = [float(m) for m in re.findall(r"cumulative = ([-\dEe.+]+)", log)]
CONT_LIMIT = 1e-4  # measured |cumulative| ~1e-7; three orders of headroom
check(
    "bounded continuity error",
    len(cum) > 0 and abs(cum[-1]) < CONT_LIMIT,
    f"final cumulative continuity = {cum[-1]:.3e} (limit {CONT_LIMIT:.0e})"
    if cum
    else "no continuity errors reported",
)

# 5. Immersed-body occupancy: lambda is a solid-volume fraction.
kind, lam = read_internal(f"{END_TIME}/lambda")
check("lambda finite", all(math.isfinite(v) for v in lam), f"{len(lam)} cells")
check(
    "lambda bounded in [0,1]",
    all(-1e-9 <= v <= 1.0 + 1e-9 for v in lam),
    f"min {min(lam):.6g}, max {max(lam):.6g}",
)
solid = [v for v in lam if v > 0.5]
# Measured: 4 cells >0.5, sum(lambda) 3.75, max 0.75. The band is wide because
# the exact count depends on how the STL cuts this coarse graded mesh; the
# point of the check is that the body is present and has not vanished or
# swallowed the domain.
check(
    "immersed body occupies cells",
    1 <= len(solid) <= 40,
    f"{len(solid)} cells with lambda>0.5 (expect 1..40, measured 4)",
)
check(
    "occupancy magnitude sane",
    1.0 <= sum(lam) <= 20.0,
    f"sum(lambda) = {sum(lam):.4g} (expect 1..20, measured 3.75)",
)

# 6. The immersed boundary is actually forcing the flow.
kind, f = read_internal(f"{END_TIME}/f")
fm = mags(f)
check("IB forcing finite", all(math.isfinite(v) for v in fm), f"{len(f)} cells")
check(
    "IB forcing non-zero",
    max(fm) > 1e-6,
    f"max|f| = {max(fm):.4g} (measured 5.08)",
)
check(
    "IB forcing bounded",
    max(fm) < 1e4,
    f"max|f| = {max(fm):.4g} < 1e4",
)

# 7. Fields are finite and have not diverged.
kind, U = read_internal(f"{END_TIME}/U")
um = mags(U)
check("U finite", all(math.isfinite(v) for v in um), f"{len(U)} cells")
# Inlet peak is 1.5; anything far above indicates divergence.
check("U bounded", max(um) < 10.0, f"max|U| = {max(um):.4g} (inlet peak 1.5)")

kind, p = read_internal(f"{END_TIME}/p")
check("p finite", all(math.isfinite(v) for v in p), f"{len(p)} cells")
check(
    "p bounded",
    max(abs(min(p)), abs(max(p))) < 1e3,
    f"p range [{min(p):.4g}, {max(p):.4g}]",
)

print()
if failures:
    print(f"FAIL: {len(failures)} of {checks} checks failed: {', '.join(failures)}")
    sys.exit(1)
print(f"PASS: all {checks} smoke checks passed")
