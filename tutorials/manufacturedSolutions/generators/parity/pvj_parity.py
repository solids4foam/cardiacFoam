"""Parity harness of the continuum junction: the C++ class pvjContinuumReference
(src/verificationModels/coupledVerification), run by the `pvjContinuumParity` OpenFOAM
application on a mesh, against the NumPy module mms_common/pvj_continuum.py.

cardiacFoam-only, like parity/test_parity.py, whose OpenFOAM environment it reuses.  The
application is built by `wmake` into parity/pvjContinuum/; `write_case` writes a polyMesh
and the case dictionary, `run_case` runs the application and parses its records
(`G` cell geometry, `S` support integrals, `I` cell integrals).
"""

import subprocess
from pathlib import Path

import numpy as np

from mms_common.pvj_meshes import polymesh_arrays
from parity.test_parity import _bash, _bashrc, parity_available  # noqa: F401  (re-exported)

HERE = Path(__file__).resolve().parent
APP_DIR = HERE / "pvjContinuum"
APP = APP_DIR / "pvjContinuumParity"
_SRC = HERE.parents[3] / "src" / "verificationModels" / "coupledVerification"
SOURCES = [
    APP_DIR / "pvjContinuumParity.C", APP_DIR / "Make" / "files", APP_DIR / "Make" / "options",
    _SRC / "pvjContinuumReference.C", _SRC / "pvjContinuumReference.H",
    HERE.parents[3] / "src" / "verificationModels" / "monodomainVerification"
    / "manufacturedFDAReference.H",
]

_HEADER = """FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
%s    object      %s;
}

"""


def build_app():
    """wmake pvjContinuumParity when missing or older than its sources."""
    if APP.exists() and all(APP.stat().st_mtime >= s.stat().st_mtime for s in SOURCES):
        return APP
    done = _bash(f'cd "{APP_DIR}" && wmake')
    if done.returncode != 0 or not APP.exists():
        raise RuntimeError(
            f"wmake of pvjContinuumParity failed\n{done.stdout[-3000:]}\n{done.stderr[-3000:]}")
    return APP


# ---- meshes: mms_common.pvj_meshes builds them; this module writes a polyMesh


def write_polymesh(case, points, cells):
    """Write constant/polyMesh of a mesh of hexahedra (cells (n, 8)) or tetrahedra (cells
    (n, 4)); one patch `walls` holds every boundary face.  Returns (faces, owner, neighbour)
    as written, for mms_common.pvj_continuum.cells_from_polymesh."""
    faces, owner, neighbour, n_internal = polymesh_arrays(points, cells)
    mesh = Path(case)/"constant"/"polyMesh"
    mesh.mkdir(parents=True, exist_ok=True)
    note = 'note        "nPoints:%d  nCells:%d  nFaces:%d  nInternalFaces:%d";\n' % (
        len(points), len(cells), len(faces), n_internal)

    def write(name, cls, body, with_note=False):
        with open(mesh/name, "w") as fh:
            fh.write(_HEADER % (cls, '    location    "constant/polyMesh";\n'
                                + ("    " + note if with_note else ""), name))
            fh.write(body)

    write("points", "vectorField", "%d\n(\n%s)\n" % (len(points), "".join(
        "(%.17g %.17g %.17g)\n" % tuple(p) for p in points)))
    write("faces", "faceList", "%d\n(\n%s)\n" % (len(faces), "".join(
        "%d(%s)\n" % (len(f), " ".join(map(str, f))) for f in faces)))
    write("owner", "labelList", "%d\n(\n%s)\n" % (len(owner), "".join("%d\n" % c for c in owner)),
          with_note=True)
    write("neighbour", "labelList", "%d\n(\n%s)\n" % (len(neighbour), "".join(
        "%d\n" % c for c in neighbour)), with_note=True)
    write("boundary", "polyBoundaryMesh",
          "1\n(\n    walls\n    {\n        type            patch;\n        nFaces          %d;\n"
          "        startFace       %d;\n    }\n)\n" % (len(faces) - n_internal, n_internal))
    return faces, owner, neighbour


def write_case(case, points, cells, kernel, radius, terminal, support, field="none", **options):
    """Write a case that pvjContinuumParity runs: the mesh and system/pvjParityDict.

    support is a dict with `kind` and, as the kind needs them, `normal` and `wallRadius`;
    options are keys read by pvjContinuumReference (referenceCutSize ...).  Returns the mesh
    as (faces, owner, neighbour)."""
    case = Path(case)
    (case/"system").mkdir(parents=True, exist_ok=True)
    mesh = write_polymesh(case, points, cells)

    def put(name, cls, body):
        with open(case/"system"/name, "w") as fh:
            fh.write(_HEADER % (cls, "", name) + body)

    put("controlDict", "dictionary",
        "application pvjContinuumParity;\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\n"
        "endTime 0;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n")
    put("fvSchemes", "dictionary",
        "ddtSchemes { default steadyState; }\ngradSchemes { default Gauss linear; }\n"
        "divSchemes { default none; }\nlaplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    put("fvSolution", "dictionary", "solvers {}\n")
    lines = ["kernel %s;" % kernel, "radius %.17g;" % radius,
             "terminal (%.17g %.17g %.17g);" % tuple(terminal), "field %s;" % field,
             "support", "{", "    kind %s;" % support["kind"]]
    if "normal" in support:
        lines.append("    normal (%.17g %.17g %.17g);" % tuple(support["normal"]))
    if "wallRadius" in support:
        lines.append("    wallRadius %.17g;" % support["wallRadius"])
    lines.append("}")
    lines += ["%s %.17g;" % (k, v) for k, v in options.items()]
    put("pvjParityDict", "dictionary", "\n".join(lines) + "\n")
    return mesh


def run_case(case):
    """Run pvjContinuumParity on a case; returns (volumes, centres, (W, F), index, omega,
    omega_f) as NumPy arrays (index: the cells that meet the ball)."""
    build_app()
    done = _bash(f'"{APP}" -case "{case}"')
    if done.returncode != 0:
        raise RuntimeError(f"pvjContinuumParity failed\n{done.stdout[-3000:]}\n{done.stderr[-3000:]}")
    geometry, support, cells = [], None, []
    for line in done.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "G" and len(parts) == 6:
            geometry.append([float(v) for v in parts[1:]])
        elif parts[0] == "S" and len(parts) == 3:
            support = (float(parts[1]), float(parts[2]))
        elif parts[0] == "I" and len(parts) == 4:
            cells.append([float(v) for v in parts[1:]])
    geometry = np.array(geometry)
    cells = np.array(cells).reshape(-1, 3)
    return (geometry[:, 1], geometry[:, 2:5], support, cells[:, 0].astype(int), cells[:, 1],
            cells[:, 2])
