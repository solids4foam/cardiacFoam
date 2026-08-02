import os
import subprocess
import re
import shutil

TPL_PATH = "three_domain_box.geo.template"
TPL_BAK = "three_domain_box.geo.template.bak"

def prepare_template(alg, netgen, smoothing):
    if not os.path.exists(TPL_BAK):
        shutil.copy(TPL_PATH, TPL_BAK)
        
    with open(TPL_BAK, "r") as f:
        base = f.read()
    
    base = re.sub(r'Mesh\.Algorithm3D = .*?;', f'Mesh.Algorithm3D = {alg};', base)
    
    additions = []
    if netgen:
        additions.append("Mesh.OptimizeNetgen = 1;")
    if smoothing:
        additions.append(f"Mesh.Smoothing = {smoothing};")
        
    base += "\n" + "\n".join(additions) + "\n"
    
    # Replace __LC__ with N=20 size (1/20 = 0.05)
    base = base.replace("__LC__", "0.05")
    
    with open("temp.geo", "w") as f:
        f.write(base)

print("Generating Original...")
prepare_template(1, 0, 0)
subprocess.run(["gmsh", "-3", "temp.geo", "-format", "vtk", "-o", "original_N20.vtk", "-v", "0"])

print("Generating Optimized...")
prepare_template(4, 1, 100)
subprocess.run(["gmsh", "-3", "temp.geo", "-format", "vtk", "-o", "optimized_N20.vtk", "-v", "0"])
