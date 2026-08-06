import os
import sys

def find_max_error():
    # Parse 0/gradError
    with open("0/gradError", "r") as f:
        lines = f.readlines()
    
    in_field = False
    in_list = False
    errors = []
    for line in lines:
        line = line.strip()
        if line.startswith("internalField"):
            if "nonuniform" in line:
                in_field = True
            elif "uniform" in line:
                val = float(line.split()[2].replace(";", ""))
                errors = [val] * 1000000 
                break
            continue
        if in_field:
            if line == "(":
                in_list = True
                continue
            if line == ")" or line == ");":
                break
            if in_list:
                try:
                    errors.append(float(line))
                except:
                    pass

    # Parse 0/C
    with open("0/C", "r") as f:
        lines = f.readlines()
    
    in_field = False
    in_list = False
    centers = []
    for line in lines:
        line = line.strip()
        if line.startswith("internalField"):
            if "nonuniform" in line:
                in_field = True
            continue
        if in_field:
            if line == "(":
                in_list = True
                continue
            if line == ")" or line == ");":
                break
            if in_list:
                try:
                    coords = line.replace("(", "").replace(")", "").split()
                    if len(coords) == 3:
                        centers.append([float(x) for x in coords])
                except:
                    pass

    max_err = -1
    max_idx = -1
    for i, err in enumerate(errors):
        if err > max_err:
            max_err = err
            max_idx = i
            
    print(f"Max Error: {max_err} at cell {max_idx}")
    if max_idx >= 0 and max_idx < len(centers):
        cx, cy, cz = centers[max_idx]
        print(f"Cell Center: x={cx:.4f}, y={cy:.4f}, z={cz:.4f}")
        # Check if it's on the boundary (domain is [0,1]^3)
        tol = 0.05
        on_bound = (cx < tol or cx > 1-tol or
                    cy < tol or cy > 1-tol or
                    cz < tol or cz > 1-tol)
        print(f"Near boundary? {on_bound}")
    else:
        print("Centers not parsed correctly or index out of range.")

if __name__ == "__main__":
    import subprocess
    N = 20
    print(f"Generating mesh for N={N}...")
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && ./Allclean > /dev/null 2>&1", shell=True, executable='/bin/bash')
    
    # Instantiate the OPTIMISED template
    lc = 1.0 / N
    with open("setup/mesh/tet/box.geo.template.optimised") as f:
        base = f.read()
    with open("box.geo", "w") as f:
        f.write(base.replace('__LC__', str(lc)))
    
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && gmsh -3 box.geo -o box.msh -format msh2 > /dev/null 2>&1", shell=True, executable='/bin/bash')
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && gmshToFoam box.msh > /dev/null 2>&1", shell=True, executable='/bin/bash')
    
    # Ensure ascii format
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && foamDictionary system/controlDict -entry writeFormat -set ascii > /dev/null 2>&1", shell=True, executable='/bin/bash')
    
    print("Running gradientReconstructionOrder...")
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && gradientReconstructionOrder > grad.log 2>&1", shell=True, executable='/bin/bash')
    print("Writing cell centers...")
    subprocess.run("source /Volumes/OpenFOAM-v2412/etc/bashrc && postProcess -func writeCellCentres > /dev/null 2>&1", shell=True, executable='/bin/bash')
    
    find_max_error()
