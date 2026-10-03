#!/usr/bin/env python3
"""Tetrahedral torso mesh around the idealized heart, conformal to it.

Run by Allmesh in electroHeartBath/meshGeneration. Reads the heart boundary
(patches EPI, BASE, ENDO_LV, ENDO_RV) from heartMesh/constant/polyMesh at
full precision, hands those triangles to gmsh
as fixed discrete surfaces, and tetrahedralises the torso box minus the heart.
Every bath triangle on the heart surface is a heart boundary face, node for
node, so stitchMesh -perfect later turns the four interface patch pairs into
internal faces.

Writes:
  torsoMesh/torso.msh         bath volume, patches torsoSurface and
                              heart_{EPI,BASE,ENDO_LV,ENDO_RV}
  constant/torsoGeometry      global-coordinate organ/box entries included
                              by system/topoSetDict

Torso frame
-----------
The box axes are the mesh axes rotated by TORSO_EULER_DEG (intrinsic x, y, z,
R = Rz Ry Rx). In that frame e1 points caudal, -e2 to the patient's left and
-e3 anterior; the heart long axis (base -> apex, mesh +x) points caudal, left
and anterior. The anterior face passes through the most anterior of
ANTERIOR_ELECTRODES and the left face through the most lateral of
LATERAL_ELECTRODES, so the precordial electrodes (the torsoECG block of
../constant/electroProperties) sit on, or at most ~9 mm beneath, the torso
surface while the heart keeps >= 15 mm from the wall.

Organs are not meshed. They are analytic shapes in the torso frame that
topoSet turns into cellZones; see system/topoSetDict.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Parameters (SI metres, degrees)
# --------------------------------------------------------------------------

TORSO_EULER_DEG = (6.0, -33.0, 18.0)

ANTERIOR_ELECTRODES = ("V1", "V2", "V3")
LATERAL_ELECTRODES = ("V4", "V5", "V6")

TORSO_DEPTH = 0.20    # anterior -> posterior (e3)
TORSO_WIDTH = 0.26    # left -> right (e2)
TORSO_LENGTH = 0.32   # cranial -> caudal (e1), centred on the heart

CHEST_WALL_THICKNESS = 0.010

# Lungs: cylinders along e1, in torso-frame coordinates (s2, s3) of the axis,
# radius, and cranial/caudal extent (s1), all relative to the heart centroid.
LUNGS = {
    "leftLung":  {"s2": -0.020, "s3": 0.100, "radius": 0.045,
                  "s1": (-0.140, 0.040)},
    "rightLung": {"s2": 0.125, "s3": 0.050, "radius": 0.060,
                  "s1": (-0.140, 0.040)},
}

# Reference point for phiE: posterior-caudal-right, inside the chest wall.
REFERENCE_INSET = 0.020

HEART_SURFACE_SIZE = 0.002   # target bath element size at the heart
TORSO_SURFACE_SIZE = 0.012   # target bath element size far from it
SIZE_GRADING_DISTANCE = 0.05

HEART_PATCHES = ("EPI", "BASE", "ENDO_LV", "ENDO_RV")


# --------------------------------------------------------------------------
# OpenFOAM ascii readers
# --------------------------------------------------------------------------

def _foam_list_body(path):
    text = Path(path).read_text()
    if "git-lfs.github.com" in text[:200]:
        sys.exit(f"ERROR: {path} is a Git-LFS pointer; run 'git lfs pull'.")
    text = text[text.index("FoamFile"):]
    text = text[text.index("}") + 1:]
    match = re.search(r"\n\s*(\d+)\s*\n\s*\(", text)
    return int(match.group(1)), text[match.end():]


def read_points(poly):
    n, body = _foam_list_body(poly / "points")
    pts = np.array(
        re.findall(r"\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)", body)[:n], dtype=float
    )
    return pts


def read_faces(poly):
    n, body = _foam_list_body(poly / "faces")
    faces = [
        [int(i) for i in m.split()]
        for m in re.findall(r"\d+\(([\d ]+)\)", body)[:n]
    ]
    return faces


def read_patches(poly):
    text = (poly / "boundary").read_text()
    patches = {}
    for name, n_faces, start in re.findall(
        r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", text
    ):
        patches[name] = (int(start), int(n_faces))
    return patches


def read_electrodes(electro_properties):
    """V-leads from the torsoECG domain's electrodePositions block."""
    text = Path(electro_properties).read_text()
    block = re.search(
        r"ecgSolver\s+torsoECG;.*?electrodePositions\s*\{(.*?)\}",
        text,
        re.S,
    )
    if block is None:
        sys.exit("ERROR: no torsoECG electrodePositions in electroProperties")
    return {
        name: np.array([float(x), float(y), float(z)])
        for name, x, y, z in re.findall(
            r"(\w+)\s*\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)", block.group(1)
        )
    }


# --------------------------------------------------------------------------
# Torso frame
# --------------------------------------------------------------------------

def rotation(euler_deg):
    a, b, c = np.radians(euler_deg)
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)],
                   [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0],
                   [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0],
                   [0, 0, 1]])
    return rz @ ry @ rx


def torso_box(heart_points, electrodes):
    """Box bounds (s_min, s_max) in torso-frame coordinates s = R^T x."""
    R = rotation(TORSO_EULER_DEG)
    heart_s = heart_points @ R
    ant = np.array([electrodes[n] for n in ANTERIOR_ELECTRODES]) @ R
    lat = np.array([electrodes[n] for n in LATERAL_ELECTRODES]) @ R

    s_min = np.empty(3)
    s_max = np.empty(3)
    s_min[2] = ant[:, 2].min()
    s_max[2] = s_min[2] + TORSO_DEPTH
    s_min[1] = lat[:, 1].min()
    s_max[1] = s_min[1] + TORSO_WIDTH
    centre1 = heart_s[:, 0].mean()
    s_min[0] = centre1 - 0.5*TORSO_LENGTH
    s_max[0] = centre1 + 0.5*TORSO_LENGTH

    all_e = np.array(list(electrodes.values())) @ R
    if np.any(all_e < s_min - 1e-12) or np.any(all_e > s_max + 1e-12):
        sys.exit("ERROR: an electrode lies outside the torso box")
    clearance = min(
        (heart_s - s_min).min(), (s_max - heart_s).min()
    )
    if clearance < CHEST_WALL_THICKNESS + 0.003:
        sys.exit(f"ERROR: heart is {clearance*1e3:.1f} mm from the torso wall")

    depth = {}
    for name, x in electrodes.items():
        s = x @ R
        depth[name] = min((s - s_min).min(), (s_max - s).min())
    return R, s_min, s_max, heart_s.mean(axis=0), clearance, depth


# --------------------------------------------------------------------------
# constant/torsoGeometry
# --------------------------------------------------------------------------

def _vec(v):
    return "({:.9g} {:.9g} {:.9g})".format(*v)


def write_torso_geometry(path, R, s_min, s_max, heart_c):
    e = R.T  # rows: e1, e2, e3 in global coordinates
    lines = [
        "// Generated by build_torso_mesh.py - torso-frame shapes in global",
        "// coordinates, included by system/topoSetDict.",
        "",
    ]

    inner_min = s_min + CHEST_WALL_THICKNESS
    inner_max = s_max - CHEST_WALL_THICKNESS
    span = inner_max - inner_min
    lines += [
        f"innerTorsoOrigin {_vec(R @ inner_min)};",
        f"innerTorsoI      {_vec(e[0]*span[0])};",
        f"innerTorsoJ      {_vec(e[1]*span[1])};",
        f"innerTorsoK      {_vec(e[2]*span[2])};",
        "",
    ]

    for name, lung in LUNGS.items():
        s2 = heart_c[1] + lung["s2"]
        s3 = heart_c[2] + lung["s3"]
        p1 = R @ np.array([heart_c[0] + lung["s1"][0], s2, s3])
        p2 = R @ np.array([heart_c[0] + lung["s1"][1], s2, s3])
        lines += [
            f"{name}P1     {_vec(p1)};",
            f"{name}P2     {_vec(p2)};",
            f"{name}Radius {lung['radius']:.9g};",
            "",
        ]

    ref = R @ (s_max - REFERENCE_INSET)
    lines += [f"phiERefPoint {_vec(ref)};", ""]
    Path(path).write_text("\n".join(lines))
    return ref


# --------------------------------------------------------------------------
# gmsh
# --------------------------------------------------------------------------

def build_msh(out_msh, points, faces, patches, R, s_min, s_max):
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("torso")

    # Box corners and faces with the built-in kernel.
    corners = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                s = np.array([
                    (s_min, s_max)[i][0],
                    (s_min, s_max)[j][1],
                    (s_min, s_max)[k][2],
                ])
                x = R @ s
                corners[i, j, k] = gmsh.model.geo.addPoint(
                    *x, TORSO_SURFACE_SIZE
                )

    def line(a, b):
        return gmsh.model.geo.addLine(corners[a], corners[b])

    edges = {}

    def edge(a, b):
        if (a, b) in edges:
            return edges[a, b]
        if (b, a) in edges:
            return -edges[b, a]
        edges[a, b] = line(a, b)
        return edges[a, b]

    box_faces = [
        [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)],
        [(1, 0, 0), (1, 0, 1), (1, 1, 1), (1, 1, 0)],
        [(0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 0, 0)],
        [(0, 1, 0), (1, 1, 0), (1, 1, 1), (0, 1, 1)],
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
        [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 1)],
    ]
    box_surfaces = []
    for quad in box_faces:
        loop = gmsh.model.geo.addCurveLoop(
            [edge(quad[n], quad[(n + 1) % 4]) for n in range(4)]
        )
        box_surfaces.append(gmsh.model.geo.addPlaneSurface([loop]))
    gmsh.model.geo.synchronize()

    # Heart boundary as discrete, pre-meshed surfaces. Node tag = OpenFOAM
    # point label + 1.
    used = sorted({p for name in HEART_PATCHES
                   for f in faces[patches[name][0]:
                                  patches[name][0] + patches[name][1]]
                   for p in f})
    heart_surfaces = {}
    for name in HEART_PATCHES:
        heart_surfaces[name] = gmsh.model.addDiscreteEntity(2)

    first = heart_surfaces[HEART_PATCHES[0]]
    gmsh.model.mesh.addNodes(
        2, first, [p + 1 for p in used], points[used].ravel().tolist()
    )
    for name in HEART_PATCHES:
        start, n = patches[name]
        tris = faces[start:start + n]
        if any(len(t) != 3 for t in tris):
            sys.exit(f"ERROR: patch {name} has non-triangular faces")
        conn = [p + 1 for t in tris for p in t]
        gmsh.model.mesh.addElementsByType(heart_surfaces[name], 2, [], conn)

    # Bath volume bounded by the box and the heart.
    volume = gmsh.model.addDiscreteEntity(
        3, -1, box_surfaces + list(heart_surfaces.values())
    )

    # Size field: fine at the heart, coarse at the torso surface.
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(
        dist, "SurfacesList", list(heart_surfaces.values())
    )
    thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thr, "InField", dist)
    gmsh.model.mesh.field.setNumber(thr, "SizeMin", HEART_SURFACE_SIZE)
    gmsh.model.mesh.field.setNumber(thr, "SizeMax", TORSO_SURFACE_SIZE)
    gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(thr, "DistMax", SIZE_GRADING_DISTANCE)
    gmsh.model.mesh.field.setAsBackgroundMesh(thr)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)

    gmsh.model.addPhysicalGroup(2, box_surfaces, name="torsoSurface")
    for name, tag in heart_surfaces.items():
        gmsh.model.addPhysicalGroup(2, [tag], name=f"heart_{name}")
    gmsh.model.addPhysicalGroup(3, [volume], name="bath")

    gmsh.model.mesh.generate(3)

    # The heart triangles must come through untouched.
    for name, tag in heart_surfaces.items():
        _, conn = gmsh.model.mesh.getElementsByType(2, tag)
        if len(conn) != 3*patches[name][1]:
            sys.exit(f"ERROR: gmsh changed the {name} triangulation")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(str(out_msh))
    n_tets = len(gmsh.model.mesh.getElementsByType(4)[0])
    gmsh.finalize()
    return n_tets


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--case-dir", default=".",
                        help="generation case (holds heartMesh/, torsoMesh/)")
    parser.add_argument("--electro-properties",
                        default="../constant/electroProperties",
                        help="electroProperties holding the torsoECG "
                             "electrodePositions, relative to --case-dir")
    parser.add_argument("--geometry-only", action="store_true",
                        help="write constant/torsoGeometry, skip gmsh")
    args = parser.parse_args()

    case = Path(args.case_dir)
    poly = case / "heartMesh" / "constant" / "polyMesh"
    points = read_points(poly)
    faces = read_faces(poly)
    patches = read_patches(poly)
    electrodes = read_electrodes(case / args.electro_properties)

    heart_ids = sorted({p for name in HEART_PATCHES
                        for f in faces[patches[name][0]:
                                       patches[name][0] + patches[name][1]]
                        for p in f})
    R, s_min, s_max, heart_c, clearance, depth = torso_box(
        points[heart_ids], electrodes
    )

    ref = write_torso_geometry(
        case / "constant" / "torsoGeometry", R, s_min, s_max, heart_c
    )

    print("Torso box (torso frame) min", s_min, "max", s_max)
    print(f"Heart-to-wall clearance {clearance*1e3:.1f} mm")
    for name, d in depth.items():
        print(f"  {name}: {d*1e3:.1f} mm beneath the torso surface")
    print("phiERefPoint", _vec(ref))

    if args.geometry_only:
        return

    n_tets = build_msh(
        case / "torsoMesh" / "torso.msh", points, faces, patches,
        R, s_min, s_max
    )
    print(f"Wrote torsoMesh/torso.msh: {n_tets} bath tetrahedra")


if __name__ == "__main__":
    main()
