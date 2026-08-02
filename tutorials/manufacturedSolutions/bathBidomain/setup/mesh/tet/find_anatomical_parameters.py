import os
import subprocess
import re
import shutil

TPL_PATH = "three_domain_box.geo.template"
TPL_BAK = "three_domain_box.geo.template.bak"

def restore_template():
    if os.path.exists(TPL_BAK):
        shutil.copy(TPL_BAK, TPL_PATH)

def run_test(alg, opt_netgen, smoothing):
    # Create new template
    with open(TPL_BAK, "r") as f:
        base = f.read()
    
    # Replace the meshing parameters at the bottom
    base = re.sub(r'Mesh\.Algorithm3D = .*?;', f'Mesh.Algorithm3D = {alg};', base)
    
    # Add smoothing and netgen if not there
    additions = []
    if opt_netgen:
        additions.append("Mesh.OptimizeNetgen = 1;")
    if smoothing > 0:
        additions.append(f"Mesh.Smoothing = {smoothing};")
        
    base += "\n" + "\n".join(additions) + "\n"
    
    with open(TPL_PATH, "w") as f:
        f.write(base)
        
    print(f"--- Testing Alg={alg}, Netgen={opt_netgen}, Smooth={smoothing} ---")
    
    # Run mesh gate N=20
    # The script uses the OpenFOAM environment, we should run it via bash
    cmd = "bash run_mesh_gate.sh 20 > /dev/null 2>&1"
    os.system(cmd)
    
    # Read checkMesh log
    log_path = "results/N20/log.checkMesh"
    if not os.path.exists(log_path):
        print("Failed to run checkMesh")
        return
        
    with open(log_path, "r") as f:
        content = f.read()
        
    m = re.search(r'Mesh non-orthogonality Max:\s+([\d\.]+)\s+average:\s+([\d\.]+)', content)
    s = re.search(r'Max skewness\s*=\s*([\d\.]+)', content)
    
    if m and s:
        max_no = float(m.group(1))
        avg_no = float(m.group(2))
        max_skew = float(s.group(1))
        print(f"Result => Mean NO: {avg_no:5.2f} | Max NO: {max_no:5.2f} | Max Skew: {max_skew:5.3f}")
    else:
        print("Could not parse checkMesh output!")

if not os.path.exists(TPL_BAK):
    shutil.copy(TPL_PATH, TPL_BAK)

try:
    # 1: Delaunay, 4: Frontal, 10: HXT
    run_test(1, 0, 0) # Baseline
    run_test(1, 1, 10)
    run_test(1, 1, 100)
    run_test(10, 0, 0)
    run_test(10, 1, 10)
    run_test(4, 1, 100)
    run_test(1, 1, 500)
finally:
    restore_template()
