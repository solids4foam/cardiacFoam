import vtk
import numpy as np
from scipy.spatial import cKDTree
import sys

def read_vtu(filepath):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.Update()
    return reader.GetOutput()

def get_point_data(mesh, array_name):
    array = mesh.GetPointData().GetArray(array_name)
    if not array:
        return None
    from vtk.util.numpy_support import vtk_to_numpy
    return vtk_to_numpy(array)

def angle_diff(v1, v2):
    # Normalize
    v1 = v1 / np.linalg.norm(v1, axis=1, keepdims=True)
    v2 = v2 / np.linalg.norm(v2, axis=1, keepdims=True)
    # Dot product (abs because fibers are bidirectional)
    dot = np.abs(np.sum(v1 * v2, axis=1))
    # Clip to valid arccos range
    dot = np.clip(dot, -1.0, 1.0)
    return np.degrees(np.arccos(dot))

def main():
    ref_path = sys.argv[1]
    of_path = sys.argv[2]
    # Scale applied to the reference points so they sit in the same physical
    # space as the OpenFOAM export. Fibers are scale-invariant; this only affects
    # the nearest-neighbour spatial match. Default 1.0 == both meshes at native
    # scale. Pass 0.02 only if comparing a native reference to a scaled export.
    ref_scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    ref_mesh = read_vtu(ref_path)
    ref_points = np.array([ref_mesh.GetPoint(i) for i in range(ref_mesh.GetNumberOfPoints())])
    ref_points_scaled = ref_points * ref_scale
    
    of_mesh = read_vtu(of_path)
    of_points = np.array([of_mesh.GetPoint(i) for i in range(of_mesh.GetNumberOfPoints())])

    print("Building KDTree mapping...")
    tree = cKDTree(ref_points_scaled)
    distances, indices = tree.query(of_points)
    print(f"Mean spatial mapping error: {np.mean(distances):.6f} m")

    ref_f0 = get_point_data(ref_mesh, 'f0')
    of_fiber = get_point_data(of_mesh, 'fiber')

    if ref_f0 is None or of_fiber is None:
        print("Error: Could not find 'f0' in reference or 'fiber' in OpenFOAM.")
        pd = of_mesh.GetPointData()
        print("OpenFOAM Point Arrays available:")
        for i in range(pd.GetNumberOfArrays()):
            print(" -", pd.GetArrayName(i))
        return

    mapped_ref_f0 = ref_f0[indices]
    
    angles = angle_diff(of_fiber, mapped_ref_f0)
    
    mean_angle = np.mean(angles)
    max_angle = np.max(angles)
    std_angle = np.std(angles)
    
    print("\n--- Fiber Validation Results ---")
    print(f"Mean Angular Error: {mean_angle:.2f} degrees")
    print(f"Max Angular Error:  {max_angle:.2f} degrees")
    print(f"Standard Dev:       {std_angle:.2f} degrees")
    
    # Save the errors back to the OpenFOAM mesh so the user can visualize where it's wrong in ParaView
    from vtk.util.numpy_support import numpy_to_vtk
    err_vtk = numpy_to_vtk(angles, deep=True)
    err_vtk.SetName("angle_error_degrees")
    of_mesh.GetPointData().AddArray(err_vtk)
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("tutorials/FiberAlignedAnisotropicMesh/simulations/monodomain_coarse_tet/VTK/error_comparison.vtu")
    writer.SetInputData(of_mesh)
    writer.Write()
    print("Wrote detailed error map to VTK/error_comparison.vtu")

if __name__ == "__main__":
    main()
