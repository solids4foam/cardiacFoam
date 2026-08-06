#!/usr/bin/env python3
"""
make_uvc_transmural.py
======================
Extracts the uvc_transmural POINT_DATA field from the VTK source file,
averages it over the 4 vertices of each tetrahedral cell, and writes the
result as an OpenFOAM volScalarField (0/uvc_transmural) in one or more
case directories.

Usage
-----
    python3 make_uvc_transmural.py

Output
------
    LBBB/0/uvc_transmural
    RBBB/0/uvc_transmural

The field is required by the ionicHeterogeneity block in electroProperties
(field name "uvc_transmural"). Values are in [0, 1]:  0 = endocardium, 1 = epicardium.
"""

import os, sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VTK_FILE   = os.path.join(SCRIPT_DIR, "LBBB", "ASCIIlegacy02_biventricular_conductivity.vtk")

CASE_DIRS  = [
    os.path.join(SCRIPT_DIR, "LBBB"),
    os.path.join(SCRIPT_DIR, "RBBB"),
]

N_CELLS    = 2_921_115
N_POINTS   = 560_974
VERTS_PER_CELL = 4         # tetrahedral mesh (offsets increase by 4)

CONN_LINE  = 511569        # first data line of CONNECTIVITY section (0-indexed)
CONN_TOTAL = N_CELLS * VERTS_PER_CELL   # 11_684_460 integers

# ── FOAMFILE header template ──────────────────────────────────────────────────
FOAM_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      uvc_transmural;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{n_cells}
(
"""

FOAM_FOOTER = """\
)
;

boundaryField
{
    defaultFaces
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""

# ── Step 1: read uvc_transmural POINT_DATA ───────────────────────────────────
print("Reading uvc_transmural POINT_DATA from VTK...")
uvc_point = np.empty(N_POINTS, dtype=np.float32)
idx = 0
with open(VTK_FILE) as f:
    in_uvc = False
    for line in f:
        if "SCALARS uvc_transmural" in line:
            in_uvc = True
            continue
        if in_uvc and line.startswith("LOOKUP_TABLE"):
            continue
        if in_uvc:
            vals = line.split()
            for v in vals:
                try:
                    uvc_point[idx] = float(v)
                    idx += 1
                    if idx == N_POINTS:
                        in_uvc = False
                        break
                except ValueError:
                    in_uvc = False
                    break
        if not in_uvc and idx == N_POINTS:
            break

print(f"  Read {idx} point values  min={uvc_point.min():.4f}  max={uvc_point.max():.4f}")

# ── Step 2: read CONNECTIVITY section ────────────────────────────────────────
print(f"Reading CONNECTIVITY ({CONN_TOTAL} integers starting at line {CONN_LINE})...")
conn_vals = np.empty(CONN_TOTAL, dtype=np.int32)
idx = 0
with open(VTK_FILE) as f:
    for lineno, line in enumerate(f):
        if lineno < CONN_LINE:
            continue
        for v in line.split():
            conn_vals[idx] = int(v)
            idx += 1
            if idx == CONN_TOTAL:
                break
        if idx == CONN_TOTAL:
            break

print(f"  Read {idx} connectivity values")
conn = conn_vals.reshape(N_CELLS, VERTS_PER_CELL)

# ── Step 3: cell-centre average ───────────────────────────────────────────────
print("Computing cell-centre averages...")
uvc_cell = uvc_point[conn].mean(axis=1)   # shape (N_CELLS,)
print(f"  Cell values  min={uvc_cell.min():.4f}  max={uvc_cell.max():.4f}  mean={uvc_cell.mean():.4f}")

# ── Step 4: write OpenFOAM field to each case ─────────────────────────────────
for case_dir in CASE_DIRS:
    out_path = os.path.join(case_dir, "0", "uvc_transmural")
    print(f"Writing {out_path} ...")
    with open(out_path, "w") as f:
        f.write(FOAM_HEADER.format(n_cells=N_CELLS))
        for v in uvc_cell:
            f.write(f"{v:.6f}\n")
        f.write(FOAM_FOOTER)
    print(f"  Done.")

print("\nAll done. Add to electroProperties monodomainSolverCoeffs:")
print("    ionicHeterogeneity")
print("    {")
print('        field             uvc_transmural;')
print("        mode              transmuralBands;")
print("        endoMInterface    0.3;")
print("        mEpiInterface     0.7;")
print("        transitionWidth   0.1;")
print("        transitionMode    blend;")
print("        smoothing         smoothstep;")
print("    }")
