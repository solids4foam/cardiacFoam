#!/usr/bin/env python3
"""Exact potentials of the lead-field ECG problems, in NumPy.

    python leadfield_reference.py --self-check
    python leadfield_reference.py --check-run RUN [--check-run RUN ...]

The slab and layered potentials repeat manufacturedLeadFieldECGReference.H
(and manufacturedFDABathPhiEAt of manufacturedFDABathBidomainReference.H with
alpha = 0 and se = 1) term by term.  The sphericalShell potential is the
transmission solution of generators/shell_mms.py, evaluated with the formula
of the zone that contains the point, as the verifier does.

--check-run reads a finished run of leadFieldECG or leadFieldECG/sphericalShell
(its log.cardiacFoam, postProcessing/manufacturedLeadFieldECG_ECG.dat and
constant/polyMesh) and compares the exactPoint and exactCentroid columns that
the verifier wrote with this module.
"""

import argparse
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATORS = os.path.abspath(os.path.join(HERE, "..", ".."))
TUTORIALS = os.path.abspath(os.path.join(GENERATORS, ".."))
SLAB_CASE = os.path.join(TUTORIALS, "leadFieldECG")
SHELL_CASE = os.path.join(SLAB_CASE, "sphericalShell")

# ------------------------------------------------------------------ slab ----

# The monodomain ECG domain of the slab case: Gi = (2 + sqrt 2) G_m and
# Ge = sqrt 2 G_m (xx entries, S/m), as constant/electroProperties
GI_XX = 0.380525375260
GE_XX = 0.157618771260
K = GI_XX/(GI_XX + GE_XX)
SIGMA_HEART = GI_XX + GE_XX

# The layered problem: interfaces, layer conductivities, surface current
X_LEFT = -0.5
X_RIGHT = 1.5
SIGMA_INNER = 0.1
SIGMA_OUTER = 0.02
ALPHA = 0.01
LAYERS = {"xLeft": X_LEFT, "xRight": X_RIGHT,
          "sigmaOuterLeft": SIGMA_OUTER, "sigmaInnerLeft": SIGMA_INNER,
          "sigmaInnerRight": SIGMA_INNER, "sigmaOuterRight": SIGMA_OUTER,
          "sigmaHeart": SIGMA_HEART}

BUILTIN_SLAB_ELECTRODES = {
    "surface1": (-0.997, 0.3137, 0.6871), "surface2": (-0.997, 0.7213, 0.2419),
    "surface3": (1.997, 0.4127, 0.5813), "surface4": (1.997, 0.8317, 0.1733),
    "surface5": (-0.3713, 0.003, 0.4419), "surface6": (-0.6871, 0.5371, 0.997),
    "surface7": (1.3317, 0.997, 0.6137), "surface8": (1.7219, 0.2713, 0.003),
    "torsoLeft": (-0.4567, 0.4137, 0.5871), "torsoLeftFar": (-0.8713, 0.2311, 0.7129),
    "torsoRight": (1.5433, 0.6173, 0.3311), "nearLeft": (-0.0213, 0.4521, 0.5317),
    "nearRight": (1.0187, 0.3719, 0.4613)}
BUILTIN_SLAB_PIN = (-0.9013, 0.0437, 0.0437)
BUILTIN_SHELL_ELECTRODES = {
    "surface1": (1.2436, 0.8151, -0.1553), "surface2": (0.2050, 0.6740, 1.3186),
    "surface3": (1.1945, -0.8631, -0.2517), "surface4": (-0.2551, 1.4080, -0.4329),
    "surface5": (-0.3258, -0.8777, 1.1656), "surface6": (0.1488, -0.3633, -1.4425),
    "surface7": (-1.4160, 0.4317, 0.2089), "surface8": (-0.8375, -1.1858, -0.3569),
    "outerNear": (0.0517, 1.0713, 0.0811), "outerSide": (-1.2117, -0.4313, 0.6271),
    "outerFar": (0.9113, 0.8121, -0.6730), "pool": (0.1131, -0.0733, 0.2117)}
BUILTIN_SHELL_PIN = (1.3137, 0.0371, 0.0459)


def slab_vm(x, t):
    """Heart Vm = sqrt(1+t) cos(pi x)."""
    return np.sqrt(1.0 + t)*np.cos(math.pi*np.asarray(x, float))


def slab_phi_raw(x, t, k=K):
    """manufacturedFDABathPhiEAt(x, t, k, 0, 1): the slab potential up to the
    gauge, with the constant C = k sqrt(1+t) of the FDA bath slab."""
    x = np.asarray(x, float)
    s = math.sqrt(1.0 + t)
    c = k*s
    return np.where(x <= 0.0, -k*s + c,
                    np.where(x <= 1.0, -k*s*np.cos(math.pi*x) + c,
                             -k*s*math.cos(math.pi) + c))


def layered_resistance(x, layers=LAYERS):
    """R(x) = integral from 0 to x of dx'/sigma(x'), branch by branch as
    manufacturedLeadFieldLayeredResistance."""
    x = np.asarray(x, float)
    xl, xr = layers["xLeft"], layers["xRight"]
    heart = 1.0/layers["sigmaHeart"]
    return np.select(
        [x <= xl, x <= 0.0, x <= 1.0, x <= xr],
        [xl/layers["sigmaInnerLeft"] + (x - xl)/layers["sigmaOuterLeft"],
         x/layers["sigmaInnerLeft"],
         x/layers["sigmaHeart"],
         heart + (x - 1.0)/layers["sigmaInnerRight"]],
        heart + (xr - 1.0)/layers["sigmaInnerRight"] + (x - xr)/layers["sigmaOuterRight"])


def layered_phi_raw(x, t, k=K, alpha=ALPHA, layers=LAYERS):
    return slab_phi_raw(x, t, k) + alpha*layered_resistance(x, layers)


def gauged(raw, x, x_ref, t, value=0.0, **kw):
    """raw(x) - raw(x_ref) + value: the gauge of a pinned conductor."""
    return raw(x, t, **kw) - raw(x_ref, t, **kw) + value


# ----------------------------------------------------------------- shell ----

_SHELL = None


def shell_model():
    """shell_mms.Model of the default parameters (built once)."""
    global _SHELL
    if _SHELL is None:
        sys.path.insert(0, GENERATORS)
        import shell_mms
        _SHELL = shell_mms.build("levelSet", shell_mms.default_params())
    return _SHELL


def shell_zone(pts):
    """0 bloodPool (|x| < r1), 1 myocardium, 2 bath (|x| > r2), as zoneOf."""
    m = shell_model()
    rho = np.sum(np.asarray(pts, float)**2, axis=-1)
    return np.where(rho < m.pv["r1"]**2, 0, np.where(rho > m.pv["r2"]**2, 2, 1))


def shell_eval(name, pts, t, zones=None, only=None):
    """`name` of shell_mms with the formula of each point's zone; points in
    zones not listed in `only` are NaN."""
    m = shell_model()
    pts = np.asarray(pts, float)
    zones = shell_zone(pts) if zones is None else np.asarray(zones)
    out = np.full(pts.shape[:-1], np.nan)
    for zone in (0, 1, 2):
        if only is not None and zone not in only:
            continue
        mask = zones == zone
        if np.any(mask):
            out[mask] = m.eval_xyz(zone, name, pts[mask], t)
    return out


def shell_phi(pts, t, x_ref, value=0.0):
    """Exact potential gauged to `value` at x_ref."""
    ref = shell_eval("phiE", np.asarray(x_ref, float)[None, :], t)[0]
    return shell_eval("phiE", pts, t) - ref + value


# ------------------------------------------------------- case dictionaries ----

def _strip(text):
    return re.sub(r"//[^\n]*", "", text)


def read_points(path, block):
    """{name: point} of the dictionary block `block` of electroProperties."""
    with open(path) as fh:
        text = _strip(fh.read())
    m = re.search(r"\b%s\s*\{([^{}]*)\}" % re.escape(block), text)
    if m is None:
        return {}
    return {g[1]: np.array([float(v) for v in g[2].split()])
            for g in re.finditer(r"(\w+)\s*\(([^)]*)\)\s*;", m.group(1))}


def read_torso_pin(path, block):
    with open(path) as fh:
        text = _strip(fh.read())
    m = re.search(r"\b%s\s*\{[^{}]*?phiERefPoint\s*\(([^)]*)\)" % re.escape(block), text)
    return None if m is None else np.array([float(v) for v in m.group(1).split()])


def case_points(kind):
    """(electrodes, pin) of the committed case; the built-in sets when absent."""
    if kind == "slab":
        path = os.path.join(SLAB_CASE, "constant", "electroProperties")
        blocks, pin_block = ("_leadFieldElectrodes",), "_leadFieldTorso"
        builtin = (BUILTIN_SLAB_ELECTRODES, BUILTIN_SLAB_PIN)
    else:
        path = os.path.join(SHELL_CASE, "constant", "electroProperties")
        blocks, pin_block = ("electrodePositions",), "torso"
        builtin = (BUILTIN_SHELL_ELECTRODES, BUILTIN_SHELL_PIN)
    if os.path.exists(path):
        electrodes = read_points(path, blocks[0])
        pin = read_torso_pin(path, pin_block)
        if electrodes and pin is not None:
            return electrodes, pin
    print("note: %s not read; using the built-in %s electrodes" % (path, kind))
    return ({k: np.array(v) for k, v in builtin[0].items()}, np.array(builtin[1]))


# ----------------------------------------------------------- self-checks ----

def check_slab_formulas():
    """Continuity of the potential and of the current, and the end values."""
    for t in (0.0, 0.05):
        s = math.sqrt(1.0 + t)
        assert abs(slab_phi_raw(-0.3, t) - 0.0) < 1e-15
        assert abs(slab_phi_raw(1.7, t) - 2*K*s) < 1e-15
        for xi in (0.0, 1.0):
            a, b = slab_phi_raw([xi - 1e-12, xi + 1e-12], t)
            assert abs(a - b) < 1e-11, ("slab continuity", xi, a, b)
        for xi in (X_LEFT, 0.0, 1.0, X_RIGHT):
            a, b = layered_resistance([xi - 1e-12, xi + 1e-12])
            assert abs(a - b) < 1e-9, ("R continuity", xi)
            # sigma R' = 1 on both sides
            h = 1e-6
            for side in (-1, 1):
                x0 = xi + side*2*h
                d = (layered_resistance(x0 + h) - layered_resistance(x0 - h))/(2*h)
                sig = {(-1, X_LEFT): SIGMA_OUTER, (1, X_LEFT): SIGMA_INNER,
                       (-1, 0.0): SIGMA_INNER, (1, 0.0): SIGMA_HEART,
                       (-1, 1.0): SIGMA_HEART, (1, 1.0): SIGMA_INNER,
                       (-1, X_RIGHT): SIGMA_INNER, (1, X_RIGHT): SIGMA_OUTER}[(side, xi)]
                assert abs(sig*d - 1.0) < 1e-8, ("sigma R'", xi, side, sig*d)
    assert abs(K - 1/math.sqrt(2)) < 1e-11


SLAB_LEVELS = (10, 20, 40, 80)
SURFACE_DEPTH = 0.003


def check_slab_electrodes():
    """The body-surface electrodes lie SURFACE_DEPTH inside the outer faces of
    the torso slabs, the probes farther inside; every electrode lies in the
    torso of both problems and off the hexahedral faces of every level."""
    e, pin = case_points("slab")
    lo, hi = np.array([-1.0, 0.0, 0.0]), np.array([2.0, 1.0, 1.0])
    for name, p in e.items():
        depth = min(np.min(p - lo), np.min(hi - p))
        if name.startswith("surface"):
            assert abs(depth - SURFACE_DEPTH) < 1e-12, (name, depth)
        else:
            assert depth > 0.01, (name, depth)
        assert p[0] < 0.0 or p[0] > 1.0, (name, "in the heart")
        for n in SLAB_LEVELS:
            # hex faces along x at -1 + i/n, along y and z at j/n
            frac = np.array([(p[0] + 1.0)*n, p[1]*n, p[2]*n])
            assert np.min(np.abs(frac - np.round(frac))) > 1e-6, (name, n)
    assert sorted(n for n in e if n.startswith("surface")) == \
        ["surface%d" % i for i in range(1, 9)]


def check_shell_zones():
    e, pin = case_points("shell")
    zones = dict(zip(e, shell_zone(np.array(list(e.values())))))
    expected = {"outerFar": 2, "outerNear": 2, "outerSide": 2, "pool": 0}
    expected.update({"surface%d" % i: 2 for i in range(1, 9)})
    assert zones == expected, zones
    assert shell_zone(pin[None, :])[0] == 2
    r = {n: float(np.linalg.norm(p)) for n, p in e.items()}
    assert all(v < 1.5 for v in r.values()), r
    assert all(abs(r["surface%d" % i] - 1.495) < 1e-4 for i in range(1, 9)), r


def run_self_checks():
    check_slab_formulas()
    check_slab_electrodes()
    check_shell_zones()


# ---------------------------------------------------------- run comparison ----

def _verifier_line(run):
    """The verifier line of log.cardiacFoam, joined with the list of electrode
    centroids that follows it when OpenFOAM prints the list over several lines."""
    with open(os.path.join(run, "log.cardiacFoam")) as fh:
        lines = fh.read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("manufacturedLeadFieldECGVerifier: reference"):
            out = [line.strip()]
            if line.rstrip().endswith("centroids"):
                for nxt in lines[i + 1:]:
                    out.append(nxt.strip())
                    if nxt.strip() == ")":
                        break
            return " ".join(out)
    raise RuntimeError("%s: no verifier line in log.cardiacFoam" % run)


def _number(line, key):
    m = re.search(r"%s ([-+0-9.eE]+)" % re.escape(key), line)
    return None if m is None else float(m.group(1))


def _samples(run):
    path = os.path.join(run, "postProcessing", "manufacturedLeadFieldECG_ECG.dat")
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.atleast_2d(np.loadtxt(path))
    return header, data


def _cells(run, points):
    """For each point: the centroid of its cell when the cell type gives it
    exactly from the vertices (tetrahedra, axis-aligned boxes), else None.
    None also where VTK places the point in no cell: VTK tests a hexahedron
    with curved (non-planar) faces by its trilinear map, OpenFOAM by its
    planar face decomposition, and the two differ next to the faceted outer
    sphere."""
    mesh = read_polymesh(run)["internalMesh"]
    xyz = read_points_file(os.path.join(run, "constant", "polyMesh", "points"))
    if len(xyz) != mesh.n_points:
        raise RuntimeError("%s: the reader added points; no exact centroids" % run)
    out = []
    for p in points:
        cell = mesh.find_containing_cell(p)
        if cell < 0:
            out.append(None)
            continue
        c = mesh.get_cell(cell)
        v = xyz[np.asarray(c.point_ids)]
        lo, hi = v.min(axis=0), v.max(axis=0)
        on_corner = (np.isclose(v, lo, rtol=0, atol=1e-12)
                     | np.isclose(v, hi, rtol=0, atol=1e-12)).all()
        exact = c.type == 10 or (c.type == 12 and bool(on_corner)
                                 and len({tuple(q) for q in np.round(v, 10)}) == 8)
        out.append(v.mean(axis=0) if exact else None)
    return out


def read_points_file(path):
    """The points of constant/polyMesh/points (ascii or binary, double)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    head = raw[:raw.index(b"}") + 1].decode("latin-1")
    body = raw[raw.index(b"}") + 1:]
    m = re.search(rb"\n\s*(\d+)\s*\n?\(", body)
    n = int(m.group(1))
    start = m.end()
    if re.search(r"format\s+binary", head):
        return np.frombuffer(body[start:start + 24*n], dtype="<f8").reshape(n, 3).copy()
    text = body[start:].decode("latin-1")
    vals = re.findall(r"\(([^()]*)\)", text)[:n]
    return np.array([[float(v) for v in q.split()] for q in vals])


def read_polymesh(path, zones=False):
    """Blocks of constant/polyMesh of the case at `path` in 64-bit floats
    ('internalMesh'; with zones, 'zones' holds the cellZones)."""
    import shutil
    import tempfile
    import pyvista as pv
    if not os.path.isdir(os.path.join(path, "constant", "polyMesh")):
        raise FileNotFoundError("%s has no constant/polyMesh" % path)
    tmp = tempfile.mkdtemp(prefix="leadfield_mesh_")
    try:
        for sub in ("constant", "system"):
            if os.path.isdir(os.path.join(path, sub)):
                os.symlink(os.path.abspath(os.path.join(path, sub)), os.path.join(tmp, sub))
        foam = os.path.join(tmp, "case.foam")
        open(foam, "w").close()
        reader = pv.OpenFOAMReader(foam)
        reader.skip_zero_time = True
        reader.reader.SetUse64BitFloats(1)
        if zones:
            reader.reader.SetReadZones(1)
        data = reader.read()
        return {name: data[name].copy() for name in data.keys()}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def check_run(run, tol=1e-12):
    """Compare the verifier columns of a run with this module; returns the
    largest differences."""
    line = _verifier_line(run)
    kind = line.split("reference ")[1].split(",")[0].strip()
    header, data = _samples(run)
    names = [h[len("exactPoint_"):] for h in header if h.startswith("exactPoint_")]
    if kind in ("slab", "layered"):
        electrodes, pin = case_points("slab")
        k = _number(line, "k")
        assert abs(k - K) < 1e-6, ("k", k, K)
        alpha = 0.0
        if kind == "layered":
            alpha = _number(line, "alpha")
            for key, ref in (("sigma outer left", SIGMA_OUTER), ("inner left", SIGMA_INNER),
                             ("heart", SIGMA_HEART), ("inner right", SIGMA_INNER),
                             ("outer right", SIGMA_OUTER)):
                assert abs(_number(line, key) - ref) <= 1e-5*ref, (key, line)

        def raw(x, t):
            x = np.asarray(x, float)[..., 0]
            return (slab_phi_raw(x, t, K) if kind == "slab"
                    else layered_phi_raw(x, t, K, alpha))
    else:
        electrodes, pin = case_points("shell")

        def raw(x, t):
            return shell_eval("phiE", np.atleast_2d(x), t)
    xe = np.array([electrodes[n] for n in names])
    centroids = _cells(run, list(xe) + [pin])
    # the log prints the centroids to 6 digits: a coarser check where the
    # vertices do not give the centroid exactly
    logged = [np.array([float(v) for v in q.split()]) for q in
              re.findall(r"\(([^()]*)\)", line.split("electrode cell centroids")[1])]
    worst = {"point": 0.0, "centroid": 0.0, "gauge": 0.0}
    logged_worst = 0.0
    for row in data:
        t = row[0]
        point = np.array([row[1 + header[1:].index("exactPoint_" + n)] for n in names])
        cent = np.array([row[1 + header[1:].index("exactCentroid_" + n)] for n in names])
        # the verifier subtracts phi at the pinned-cell centroid: one constant
        gauge = raw(xe, t) - point
        worst["point"] = max(worst["point"], float(np.ptp(gauge)))
        if centroids[-1] is not None:
            worst["gauge"] = max(worst["gauge"],
                                 float(abs(gauge[0] - raw(centroids[-1][None, :], t)[0])))
        for i, c in enumerate(centroids[:-1]):
            if c is not None:
                worst["centroid"] = max(worst["centroid"],
                                        float(abs(raw(c[None, :], t)[0] - gauge[0] - cent[i])))
            else:
                logged_worst = max(logged_worst, float(abs(
                    raw(logged[i][None, :], t)[0] - gauge[0] - cent[i])))
    exact_cells = sum(c is not None for c in centroids)
    print("%s: reference %s, %d samples, %d electrodes; spread of the gauge over the "
          "electrodes %.2e; gauge against phi at the pinned-cell centroid %s; "
          "exactCentroid %s%s"
          % (run, kind, len(data), len(names), worst["point"],
             "%.2e" % worst["gauge"] if centroids[-1] is not None else "not checked",
             "%.2e (%d of %d cells)" % (worst["centroid"], exact_cells - (centroids[-1] is not None),
                                         len(names)),
             "; exactCentroid at the logged 6-digit centroids %.2e" % logged_worst
             if exact_cells < len(centroids) else ""))
    assert all(v < tol for v in worst.values()), (run, worst)
    assert logged_worst < 1e-4, (run, logged_worst)
    return worst


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--self-check", action="store_true")
    ap.add_argument("--check-run", action="append", default=[], metavar="RUN")
    args = ap.parse_args(argv)
    if not args.self_check and not args.check_run:
        ap.error("give --self-check or --check-run RUN")
    if args.self_check:
        run_self_checks()
        print("self-checks passed")
    for run in args.check_run:
        check_run(run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
