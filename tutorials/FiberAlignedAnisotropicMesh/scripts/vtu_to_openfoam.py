import sys
import os
import re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.spatial import cKDTree

def parse_openfoam_vector_field(filepath):
    """Parses an OpenFOAM volVectorField and returns a numpy array of the internalField."""
    with open(filepath, 'r') as f:
        content = f.read()

    # Find the internalField nonuniform List<vector>
    match = re.search(r'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;', content, re.DOTALL)
    if not match:
        raise ValueError(f"Could not parse internalField from {filepath}")
    
    num_elements = int(match.group(1))
    data_str = match.group(2)
    
    # Extract all numbers
    numbers = np.fromiter((float(x) for x in re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+', data_str)), dtype=np.float64)
    vectors = numbers.reshape((num_elements, 3))
    return vectors

def write_openfoam_vector_field(out_path, field_name, vectors, boundary_patches):
    """Writes a volVectorField to OpenFOAM format."""
    num_elements = vectors.shape[0]
    
    with open(out_path, 'w') as f:
        f.write("FoamFile\n{\n")
        f.write("    version     2.0;\n")
        f.write("    format      ascii;\n")
        f.write("    class       volVectorField;\n")
        f.write('    location    "0";\n')
        f.write(f'    object      {field_name};\n')
        f.write("}\n")
        f.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")
        f.write("dimensions      [0 0 0 0 0 0 0];\n\n")
        f.write(f"internalField   nonuniform List<vector> \n{num_elements}\n(\n")
        
        # Write vectors
        for i in range(num_elements):
            f.write(f"({vectors[i,0]:.6g} {vectors[i,1]:.6g} {vectors[i,2]:.6g})\n")
            
        f.write(")\n;\n\n")
        
        f.write("boundaryField\n{\n")
        for patch in boundary_patches:
            f.write(f"    {patch}\n    {{\n")
            f.write("        type            calculated;\n")
            f.write("        value           uniform (0 0 0);\n")
            f.write("    }\n")
        f.write("}\n")
        f.write("// ************************************************************************* //\n")

def get_vtk_point_data(mesh, name):
    arr = mesh.GetPointData().GetArray(name)
    if arr:
        return vtk_to_numpy(arr)
    return None

def main():
    if len(sys.argv) < 4:
        print("Usage: python vtu_to_openfoam.py <matlab_result.vtu> <openfoam_0_dir> <scale_factor>")
        sys.exit(1)
        
    vtu_path = sys.argv[1]
    of_0_dir = sys.argv[2]
    scale_factor = float(sys.argv[3]) # Scale for VTU to match OpenFOAM
    
    # 1. Read OpenFOAM Cell Centers
    c_path = os.path.join(of_0_dir, 'C')
    if not os.path.exists(c_path):
        print(f"Error: {c_path} does not exist. Run 'writeCellCentres' in OpenFOAM first.")
        sys.exit(1)
        
    print(f"Parsing OpenFOAM cell centers from {c_path}...")
    of_centers = parse_openfoam_vector_field(c_path)
    print(f"Loaded {of_centers.shape[0]} cell centers from OpenFOAM.")
    
    # 2. Extract Boundary Patches from C file to reuse them
    with open(c_path, 'r') as f:
        c_content = f.read()
    
    boundary_patches = []
    # Find everything inside boundaryField { ... }
    bf_match = re.search(r'boundaryField\s*\{([^}]*)\}', c_content)
    if bf_match:
        bf_block = bf_match.group(1)
        # Find all top-level patch names (words before '{')
        # This is a simple regex that usually works for simple boundary files
        patches = re.findall(r'^\s*([A-Za-z0-9_]+)\s*\{', bf_block, re.MULTILINE)
        boundary_patches = patches
    
    print(f"Found boundary patches: {boundary_patches}")
    
    # 3. Load MATLAB VTU
    print(f"Loading MATLAB VTU: {vtu_path}...")
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    vtu_mesh = reader.GetOutput()
    
    # Get VTU points and scale them
    vtu_points = np.array([vtu_mesh.GetPoint(i) for i in range(vtu_mesh.GetNumberOfPoints())])
    vtu_points_scaled = vtu_points * scale_factor
    
    # Get fields from MATLAB
    f0 = get_vtk_point_data(vtu_mesh, 'f0')
    s0 = get_vtk_point_data(vtu_mesh, 's0')
    n0 = get_vtk_point_data(vtu_mesh, 'n0')
    
    if f0 is None:
        print("Error: 'f0' not found in MATLAB VTU.")
        sys.exit(1)
        
    print("Building KDTree of VTU points...")
    tree = cKDTree(vtu_points_scaled)
    distances, indices = tree.query(of_centers)
    print(f"Mean spatial mapping error: {np.mean(distances):.6e} m")
    
    print("Mapping fields...")
    of_fiber = f0[indices]
    
    # Write fiber
    out_fiber = os.path.join(of_0_dir, 'fiber')
    write_openfoam_vector_field(out_fiber, 'fiber', of_fiber, boundary_patches)
    print(f"Wrote {out_fiber}")
    
    if s0 is not None:
        of_sheet = s0[indices]
        out_sheet = os.path.join(of_0_dir, 'sheet')
        write_openfoam_vector_field(out_sheet, 'sheet', of_sheet, boundary_patches)
        print(f"Wrote {out_sheet}")
        
    if n0 is not None:
        of_normal = n0[indices]
        out_normal = os.path.join(of_0_dir, 'normalDirection')
        write_openfoam_vector_field(out_normal, 'normalDirection', of_normal, boundary_patches)
        print(f"Wrote {out_normal}")
        
if __name__ == "__main__":
    main()
