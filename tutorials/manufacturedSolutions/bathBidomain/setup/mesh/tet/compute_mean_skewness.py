#!/usr/bin/env python3
import os
import re
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CASE_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../../.."))
TPL_OPT = os.path.join(SCRIPT_DIR, "three_domain_box.geo.template.optimised")
GEO_WORK = os.path.join(SCRIPT_DIR, "three_domain_box.geo")

OF_BASHRC = "/Volumes/OpenFOAM-v2412/etc/bashrc"
env_cmd = f"source {OF_BASHRC} > /dev/null 2>&1 && "

def shell(cmd, cwd=CASE_ROOT):
    return subprocess.run(
        ["bash", "-c", env_cmd + cmd],
        capture_output=True, text=True, cwd=cwd
    )

def read_vol_scalar_field(path):
    if not os.path.exists(path): return None
    with open(path, "r", errors="replace") as fh: text = fh.read()
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", text)
    if m: return [float(m.group(1))]
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\(", text)
    if not m: return None
    n = int(m.group(1))
    start = text.index("(", m.end() - 1) + 1
    end = text.index(")", start)
    vals = text[start:end].split()
    return [float(v) for v in vals]

with open(TPL_OPT) as f:
    base = f.read()

for N in [10, 20, 40, 80]:
    lc = 1.0 / N
    geo_text = re.sub(r'__LC__', str(lc), base)
    with open(GEO_WORK, "w") as f:
        f.write(geo_text)
    
    msh = os.path.join(CASE_ROOT, "three_domain_box.msh")
    print(f"Generating N={N}...")
    shell(f"gmsh -3 {GEO_WORK} -o {msh} -format msh2 -v 0")
    shell(f"gmshToFoam {msh}")
    os.remove(msh)
    
    shell("checkMesh -allGeometry -allTopology -writeAllFields -time 0")
    
    skew_path = os.path.join(CASE_ROOT, "0", "skewness")
    skew_vals = read_vol_scalar_field(skew_path)
    if skew_vals:
        mean_skew = sum(skew_vals) / len(skew_vals)
        max_skew = max(skew_vals)
        print(f"N={N}: Mean skewness = {mean_skew:.6f}, Max skewness = {max_skew:.6f}")
    else:
        print(f"N={N}: Failed to read 0/skewness")
    
    shell("./Allclean")
