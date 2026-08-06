import os
import subprocess
import re
import shutil

TPL_PATH = "three_domain_box.geo.template"
TPL_BAK = "three_domain_box.geo.template.bak"

def prepare_template():
    if not os.path.exists(TPL_BAK):
        shutil.copy(TPL_PATH, TPL_BAK)
        
    with open(TPL_BAK, "r") as f:
        base = f.read()
    
    # Replace the meshing parameters
    base = re.sub(r'Mesh\.Algorithm3D = .*?;', 'Mesh.Algorithm3D = 4;', base)
    
    additions = [
        "Mesh.OptimizeNetgen = 1;",
        "Mesh.Smoothing = 100;"
    ]
    base += "\n" + "\n".join(additions) + "\n"
    
    with open(TPL_PATH, "w") as f:
        f.write(base)

def run_ladder():
    prepare_template()
    
    results = {}
    try:
        for N in [10, 20, 40, 80]:
            print(f"--- Running N={N} ---")
            os.system(f"bash run_mesh_gate.sh {N} > /dev/null 2>&1")
            
            log_path = f"results/N{N}/log.checkMesh"
            if not os.path.exists(log_path):
                print(f"Failed to run checkMesh for N={N}")
                continue
                
            with open(log_path, "r") as f:
                content = f.read()
                
            m = re.search(r'Mesh non-orthogonality Max:\s+([\d\.]+)\s+average:\s+([\d\.]+)', content)
            s = re.search(r'Max skewness\s*=\s*([\d\.]+)', content)
            
            if m and s:
                max_no = float(m.group(1))
                avg_no = float(m.group(2))
                max_skew = float(s.group(1))
                results[N] = (avg_no, max_no, max_skew)
                print(f"N={N:2d} => Mean NO: {avg_no:5.2f} | Max NO: {max_no:5.2f} | Max Skew: {max_skew:5.3f}")
            else:
                print(f"Could not parse N={N}")
    finally:
        if os.path.exists(TPL_BAK):
            shutil.copy(TPL_BAK, TPL_PATH)
            
    # Write markdown table summary
    with open("ladder_summary.md", "w") as f:
        f.write("| Resolution (N) | Mean NO | Max NO | Max Skew |\n")
        f.write("|----------------|---------|--------|----------|\n")
        for N in [10, 20, 40, 80]:
            if N in results:
                avg, mx, sk = results[N]
                f.write(f"| {N:14d} | {avg:7.2f}° | {mx:6.2f}° | {sk:8.3f} |\n")

run_ladder()
