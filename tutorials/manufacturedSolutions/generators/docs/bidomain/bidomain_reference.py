#!/usr/bin/env python3
"""Bidomain manufactured solutions: case checks, run comparisons and figure helpers.

    python bidomain_reference.py --self-check
    python bidomain_reference.py --run-check fda=CASE ... shell=CASE ...

The exact fields come from the generators: reference/fda_bidomain.py, emitted
by fda_bidomain_mms.py, for the FDA cube, and shell_bidomain_mms.py for the
sealed shell.  This module adds what depends on the tutorial cases.

--self-check asserts that the conductivities and phiERefPoint of the
tutorials manufacturedSolutions/bidomain and bidomain/sphericalShell equal
those of the generators.  --run-check reads the cell centres (computed from
constant/polyMesh as OpenFOAM does), the fields of the last time directory
and postProcessing/*.dat of finished runs and checks that the exact fields
reproduce the errors the cardiacFoam verifier reported; it writes nothing.
"""

import argparse
import glob
import importlib.util
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GENERATORS = os.path.abspath(os.path.join(HERE, "..", ".."))
TUTORIAL = os.path.abspath(os.path.join(GENERATORS, "..", "bidomain"))
sys.path.insert(0, GENERATORS)

import shell_bidomain_mms as sbm  # noqa: E402


def _load(name):
    path = os.path.join(GENERATORS, "reference", name + ".py")
    spec = importlib.util.spec_from_file_location("ref_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FDA = _load("fda_bidomain")

# ------------------------------------------------------------------ FDA cube ----

GI = tuple(FDA.PARAMS["GI_DIAG"])
GE = tuple(FDA.PARAMS["GE_DIAG"])
PHIE_REF_POINT = tuple(FDA.PARAMS["PHIE_REF_POINT"])
PHIE_REF_VALUE = FDA.PARAMS["PHIE_REF_VALUE"]


def kappa():
    return FDA.PARAMS["KAPPA"]


def beta(d):
    return FDA.PARAMS["BETA"][d]


def _field(name, d, x, y, z, t):
    return getattr(FDA, "%s_d%d" % (name, d))(x, y, z, t)


def Vm(d, x, y, z, t):
    return _field("Vm", d, x, y, z, t)


def phiE(d, x, y, z, t):
    return _field("phiE", d, x, y, z, t)


def phiI(d, x, y, z, t):
    return _field("phiI", d, x, y, z, t)


def u1(d, x, y, z, t):
    return _field("u1", d, x, y, z, t)


def read_dictionary_conductivities(path):
    text = open(path).read()
    out = {}
    for name in ("conductivityIntracellular", "conductivityExtracellular"):
        m = re.search(name + r"\s*\[[^\]]*\]\s*\(([^)]*)\)", text)
        out[name] = [float(v) for v in m.group(1).split()]
    m = re.search(r"phiERefPoint\s*\(([^)]*)\)", text)
    out["phiERefPoint"] = tuple(float(v) for v in m.group(1).split())
    return out


def fda_case_checks():
    case = read_dictionary_conductivities(os.path.join(TUTORIAL, "constant", "electroProperties"))
    gi6, ge6 = case["conductivityIntracellular"], case["conductivityExtracellular"]
    assert gi6 == [GI[0], 0, 0, GI[1], 0, GI[2]], "Gi of the tutorial"
    assert ge6 == [GE[0], 0, 0, GE[1], 0, GE[2]], "Ge of the tutorial"
    assert case["phiERefPoint"] == PHIE_REF_POINT, "phiERefPoint of the tutorial"


# ------------------------------------------------------------- sealed shell ----

class Shell:
    """Adapter with the interface of the emitted shell reference modules
    (PARAMS, field(name, x, y, z)) used by docs/_common/figstyle.py."""

    PARAMS = {"default": sbm.default_params()}

    @staticmethod
    def evaluate(name, x, y, z):
        pts = np.stack(np.broadcast_arrays(*(np.asarray(v, float) for v in (x, y, z))), -1)
        return sbm.evaluate(name, pts)

    def field(self, name, x, y, z):
        if name == "sigma":
            return self.tensor("Gi", x, y, z)
        if name == "Vm":
            return self.evaluate("Fi", x, y, z) - self.evaluate("Fe", x, y, z)
        return self.evaluate(name, x, y, z)

    def tensor(self, name, x, y, z):
        """(3, 3, N) array of Gi or Ge."""
        T = sbm._full(self.evaluate(name, x, y, z))
        return np.moveaxis(T.reshape(-1, 3, 3), 0, -1)


SHELL = Shell()


def shell_volume_mean(name, nr=40, nc=40, nphi=80):
    """Volume mean over the shell of a scalar of shell_bidomain_mms."""
    p = sbm.default_params()
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xc, wc = np.polynomial.legendre.leggauss(nc)
    r = 0.5*(p["r2"] - p["r1"])*xr + 0.5*(p["r2"] + p["r1"])
    wr = 0.5*(p["r2"] - p["r1"])*wr*r**2
    ph = 2*math.pi*(np.arange(nphi) + 0.5)/nphi
    R, C, PH = np.meshgrid(r, xc, ph, indexing="ij")
    S = np.sqrt(1 - C**2)
    pts = np.stack([R*S*np.cos(PH), R*S*np.sin(PH), R*C], -1)
    W = np.broadcast_to(wr[:, None, None]*wc[None, :, None]*(2*math.pi/nphi), R.shape)
    volume = 4.0*math.pi*(p["r2"]**3 - p["r1"]**3)/3.0
    assert abs(W.sum() - volume) < 1e-12*volume, "shell quadrature weights"
    return float((W*sbm.evaluate(name, pts)).sum()/W.sum())


def read_shell_dictionary():
    text = open(os.path.join(TUTORIAL, "sphericalShell", "constant", "electroProperties")).read()
    m = re.search(r"phiERefPoint\s*\(([^)]*)\)", text)
    return {"phiERefPoint": tuple(float(v) for v in m.group(1).split()),
            "ionicModel": re.search(r"ionicModel\s+(\w+)", text).group(1)}


def shell_case_checks():
    case = read_shell_dictionary()
    assert case["phiERefPoint"] == sbm.PHIE_REF_POINT, "phiERefPoint of the shell tutorial"
    assert case["ionicModel"] == "monodomainFDAManufactured", "ionic model of the shell tutorial"


def run_self_checks():
    fda_case_checks()
    shell_case_checks()
    print("the tutorial dictionaries match the generators")


# ------------------------------------------------------- comparison with runs ----

def _foam_lists(path):
    """Every list (count followed by '(') of an OpenFOAM file after its header,
    as raw bytes, in order; binary or ascii."""
    data = open(path, "rb").read()
    fmt = re.search(rb"format\s+(\w+);", data).group(1).decode()
    cls = re.search(rb"class\s+(\w+);", data).group(1).decode()
    body = data[data.index(b"}", data.index(b"FoamFile")) + 1:]
    return fmt, cls, body


_WIDTH = {"scalar": 1, "vector": 3, "symmTensor": 6, "tensor": 9}


def read_internal_field(path):
    fmt, cls, body = _foam_lists(path)
    kind = cls.replace("vol", "").replace("Field", "")
    kind = kind[0].lower() + kind[1:]
    m = re.search(rb"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", body)
    if m is None:
        u = re.search(rb"internalField\s+uniform\s+([^;]+);", body)
        return float(u.group(1))
    n, w = int(m.group(1)), _WIDTH[kind]
    start = m.end()
    if fmt == "binary":
        a = np.frombuffer(body[start:start + 8*n*w], dtype="<f8")
    else:
        a = np.array(body[start:body.index(b")\n;", start)].replace(b"(", b" ").replace(b")", b" ")
                     .split(), float)
    return a.reshape(n, w) if w > 1 else a


def _read_list(path, dtype, width=1, which=0):
    fmt, cls, body = _foam_lists(path)
    pos, k = 0, 0
    while True:
        m = re.compile(rb"(\d+)\s*\(").search(body, pos)
        n = int(m.group(1))
        start = m.end()
        size = n*width*np.dtype(dtype).itemsize
        if fmt == "binary":
            a = np.frombuffer(body[start:start + size], dtype=dtype)
            end = start + size
        else:
            end = body.index(b"\n)", start)
            text = body[start:end].replace(b"(", b" ").replace(b")", b" ").split()
            a = np.array(text, dtype=dtype)
        if k == which:
            return a.reshape(n, width) if width > 1 else a
        pos, k = end + 1, k + 1


def cell_centres_and_volumes(case):
    """Cell centres and volumes of constant/polyMesh, by the face-triangle and
    face-pyramid decomposition of OpenFOAM's primitiveMeshTools."""
    pm = os.path.join(case, "constant", "polyMesh")
    p = _read_list(os.path.join(pm, "points"), "<f8", 3)
    fmt, cls, _ = _foam_lists(os.path.join(pm, "faces"))
    if cls == "faceCompactList":
        offs = _read_list(os.path.join(pm, "faces"), "<i4", 1, 0)
        labs = _read_list(os.path.join(pm, "faces"), "<i4", 1, 1)
        faces = [labs[offs[i]:offs[i + 1]] for i in range(len(offs) - 1)]
    else:
        raise NotImplementedError("faceList")
    own = _read_list(os.path.join(pm, "owner"), "<i4")
    nei = _read_list(os.path.join(pm, "neighbour"), "<i4")
    nf = len(faces)
    fc = np.zeros((nf, 3))
    fa = np.zeros((nf, 3))
    for i, f in enumerate(faces):
        q = p[f]
        if len(f) == 3:
            fc[i] = q.mean(0)
            fa[i] = 0.5*np.cross(q[1] - q[0], q[2] - q[0])
            continue
        c0 = q.mean(0)
        nxt = np.roll(q, -1, axis=0)
        n = np.cross(nxt - q, c0 - q)
        a = np.linalg.norm(n, axis=1)
        fc[i] = (1.0/3.0)*((a[:, None]*(q + nxt + c0)).sum(0))/a.sum()
        fa[i] = 0.5*n.sum(0)
    ncell = int(max(own.max(), nei.max() if len(nei) else 0)) + 1
    cest = np.zeros((ncell, 3))
    cnt = np.zeros(ncell)
    np.add.at(cest, own, fc)
    np.add.at(cnt, own, 1)
    np.add.at(cest, nei, fc[:len(nei)])
    np.add.at(cnt, nei, 1)
    cest /= cnt[:, None]
    cc = np.zeros((ncell, 3))
    vol = np.zeros(ncell)
    pyr = np.einsum("ij,ij->i", fa, fc - cest[own])
    np.add.at(cc, own, pyr[:, None]*(0.75*fc + 0.25*cest[own]))
    np.add.at(vol, own, pyr)
    ni = np.arange(len(nei))
    pyr = np.einsum("ij,ij->i", fa[ni], cest[nei] - fc[ni])
    np.add.at(cc, nei, pyr[:, None]*(0.75*fc[ni] + 0.25*cest[nei]))
    np.add.at(vol, nei, pyr)
    cc /= vol[:, None]
    return cc, vol/3.0


def _norms(e, V):
    e = np.abs(e)
    return (float((V*e).sum()/V.sum()), float(np.sqrt((V*e*e).sum()/V.sum())), float(e.max()))


def _rows(dat):
    out = {}
    for line in open(dat):
        parts = line.split()
        if len(parts) == 4 and not line.startswith("#"):
            try:
                out[parts[0]] = tuple(float(v) for v in parts[1:])
            except ValueError:
                pass
    return out


def _last_time(case):
    times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and d != "0"]
    return max(times, key=float)


def _compare(label, mine, reported, rtol):
    worst = max(abs(a - b)/max(abs(b), 1e-300) for a, b in zip(mine, reported))
    ok = worst <= rtol
    print("  %-12s computed %s reported %s  rel. diff %.1e %s" % (
        label, " ".join("%.6g" % v for v in mine), " ".join("%.6g" % v for v in reported),
        worst, "ok" if ok else "FAIL"))
    return ok


def run_check_fda(case):
    """Vm, phiE (with the reference-point constant), phiE_gauge, phiI_gauge
    errors and the gauge offset of an FDA run, recomputed from its fields."""
    T = _last_time(case)
    t = float(T)
    dat = glob.glob(os.path.join(case, "postProcessing", "*_cells.dat"))[0]
    d = int(os.path.basename(dat)[0])
    rows = _rows(dat)
    offset = float(re.search(r"phiE gauge offset\s*=\s*(\S+)", open(dat).read()).group(1))
    C, V = cell_centres_and_volumes(case)
    x, y, z = C.T
    vm = read_internal_field(os.path.join(case, T, "Vm"))
    pe = read_internal_field(os.path.join(case, T, "phiE"))
    u1_0 = read_internal_field(os.path.join(case, T, "preProcess_u1"))
    ex_v, ex_e = Vm(d, x, y, z, t), phiE(d, x, y, z, t)
    off = float((V*(pe - ex_e)).sum()/V.sum())
    print("FDA %s (d=%d, t=%g, %d cells)" % (os.path.basename(case.rstrip("/")), d, t, len(V)))
    ok = _compare("Vm", _norms(vm - ex_v, V), rows["Vm"], 2e-5)
    ok &= _compare("phiE", _norms(pe - ex_e, V), rows["phiE"], 2e-5)
    ok &= _compare("gauge offset", (off,), (offset,), 2e-5)
    ok &= _compare("phiE_gauge", _norms(pe - ex_e - off, V), rows["phiE_gauge"], 2e-5)
    ok &= _compare("phiI_gauge", _norms(vm + pe - (ex_v + ex_e) - off, V), rows["phiI_gauge"], 2e-5)
    e0 = float(np.abs(u1_0 - u1(d, x, y, z, 0.0)).max())
    print("  preProcess_u1 against u1(t=0): max |diff| %.1e" % e0)
    return ok and e0 < 1e-12


def run_check_shell(case):
    """Gi, Ge written by setManufacturedFields, and the Vm, phiE_gauge errors
    of a shell run, recomputed from its fields."""
    T = _last_time(case)
    t = float(T)
    dat = glob.glob(os.path.join(case, "postProcessing", "*_cells.dat"))[0]
    rows = _rows(dat)
    C, V = cell_centres_and_volumes(case)
    print("shell %s (t=%g, %d cells)" % (os.path.basename(case.rstrip("/")), t, len(V)))
    ok = True
    for name, fname in (("Gi", "ConductivityIntracellular"), ("Ge", "ConductivityExtracellular")):
        G = read_internal_field(os.path.join(case, "0", fname))
        ref = sbm.evaluate(name, C)
        rel = float((np.linalg.norm(G - ref, axis=1)/np.linalg.norm(ref, axis=1)).max())
        print("  0/%s against %s at the cell centres: max rel. diff %.1e" % (fname, name, rel))
        ok &= rel < 1e-9
    S = math.sqrt(1 + t)
    Fi, Fe = sbm.evaluate("Fi", C), sbm.evaluate("Fe", C)
    ref = np.array(read_shell_dictionary()["phiERefPoint"])
    shift = -S*float(sbm.evaluate("Fe", ref))
    vm = read_internal_field(os.path.join(case, T, "Vm"))
    pe = read_internal_field(os.path.join(case, T, "phiE"))
    ex_e = S*Fe + shift
    off = float((V*(pe - ex_e)).sum()/V.sum())
    ok &= _compare("Vm", _norms(vm - S*(Fi - Fe), V), rows["Vm"], 1e-8)
    ok &= _compare("phiE_gauge", _norms(pe - ex_e - off, V), rows["phiE_gauge"], 1e-8)
    u1_0 = read_internal_field(os.path.join(case, T, "preProcess_u1"))
    g = sbm.evaluate("g", C)
    e0 = float(np.abs(u1_0 - (g + (Fi - Fe))).max())
    print("  preProcess_u1 against u1(t=0): max |diff| %.1e" % e0)
    return ok and e0 < 1e-12


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--self-check", action="store_true")
    grp.add_argument("--run-check", nargs="+", metavar="KIND=CASE",
                     help="KIND is fda or shell; CASE a finished run")
    args = ap.parse_args(argv)
    if args.self_check:
        run_self_checks()
        return 0
    ok = True
    for item in args.run_check:
        kind, case = item.split("=", 1)
        ok &= {"fda": run_check_fda, "shell": run_check_shell}[kind](case)
    print("run checks passed" if ok else "run checks FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
