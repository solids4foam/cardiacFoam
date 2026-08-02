#!/usr/bin/env python3
"""
Characterize mesh quality of the anatomically-optimised single-box tet family
for the eikonal MMS case, before committing to full solver runs.

Runs gmsh + gmshToFoam + checkMesh from the eikonal case root (where
system/controlDict etc. live), exactly as run_eikonal_tet.sh does.
"""
import os
import re
import subprocess
import sys

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
CASE_ROOT   = os.path.abspath(os.path.join(SCRIPT_DIR, "../../.."))
TPL_ORIG    = os.path.join(SCRIPT_DIR, "box.geo.template")
TPL_OPT     = os.path.join(SCRIPT_DIR, "box.geo.template.optimised")
GEO_WORK    = os.path.join(SCRIPT_DIR, "box.geo")   # where run scripts put it

# Source OpenFOAM
OF_BASHRC = "/Volumes/OpenFOAM-v2412/etc/bashrc"
env_cmd   = f"source {OF_BASHRC} > /dev/null 2>&1 && "

def shell(cmd, cwd=CASE_ROOT):
    return subprocess.run(
        ["bash", "-c", env_cmd + cmd],
        capture_output=True, text=True, cwd=cwd
    )

# ── Build optimised template once ────────────────────────────────────────────
with open(TPL_ORIG) as f:
    base = f.read()

base = re.sub(
    r'Mesh\.Algorithm3D\s*=\s*\d+\s*;.*',
    'Mesh.Algorithm3D = 4;  // Frontal: anatomically-matched quality',
    base
)
base = base.rstrip() + "\nMesh.OptimizeNetgen = 1;\nMesh.Smoothing = 100;\n"

with open(TPL_OPT, "w") as f:
    f.write(base)

print(f"Optimised template written to:\n  {TPL_OPT}\n")
print(f"{'N':>4}  {'Mean NO':>8}  {'Max NO':>8}  {'Max Skew':>10}")
print("-" * 40)

results = {}
for N in [10, 20, 40, 80]:
    lc = 1.0 / N
    geo_text = re.sub(r'__LC__', str(lc), base)

    # Write geo to the location the run scripts use
    with open(GEO_WORK, "w") as f:
        f.write(geo_text)

    msh = os.path.join(CASE_ROOT, "box.msh")

    # 1. gmsh
    r = shell(f"gmsh -3 {GEO_WORK} -o {msh} -format msh2 -v 0")
    if r.returncode != 0 or not os.path.exists(msh):
        print(f"  N={N}: gmsh FAILED"); print(r.stderr[-300:]); continue

    # 2. gmshToFoam  (must run from case root — that's where system/ lives)
    r = shell(f"gmshToFoam {msh}")
    os.remove(msh)
    if r.returncode != 0:
        print(f"  N={N}: gmshToFoam FAILED"); print(r.stderr[-300:]); continue

    # 3. checkMesh
    r = shell("checkMesh -allGeometry -allTopology")
    log = r.stdout + r.stderr

    m = re.search(r'Mesh non-orthogonality Max:\s+([\d.]+)\s+average:\s+([\d.]+)', log)
    s = re.search(r'Max skewness\s*=\s*([\d.]+)', log)

    if m and s:
        max_no  = float(m.group(1))
        mean_no = float(m.group(2))
        max_sk  = float(s.group(1))
        results[N] = (mean_no, max_no, max_sk)
        print(f"  N={N:2d}  {mean_no:8.2f}°  {max_no:8.2f}°  {max_sk:10.3f}")
    else:
        print(f"  N={N}: could not parse checkMesh output")
        for line in log.splitlines()[-20:]:
            print("   |", line)

    # Clean up polyMesh so next N starts fresh
    shell("./Allclean")

print("-" * 40)
print("Target: Mean NO ~14-16° (anatomical range)")
