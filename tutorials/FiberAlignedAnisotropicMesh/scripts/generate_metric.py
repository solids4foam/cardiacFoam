import sys
import numpy as np
import vtk
from vtk.util import numpy_support
import meshio
import argparse

def generate_metric(input_vtu, output_msh, h_f=0.6, h_t=0.15, alpha_endo=60.0, alpha_epi=-60.0):
    print(f"Reading mesh from {input_vtu}...")
    
    # Read the VTU file using VTK to compute gradients
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_vtu)
    reader.Update()
    vtk_mesh = reader.GetOutput()
    
    # Compute gradient of 'tm' (transmural)
    grad_filter = vtk.vtkGradientFilter()
    grad_filter.SetInputData(vtk_mesh)
    grad_filter.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "tm")
    grad_filter.SetResultArrayName("grad_tm")
    grad_filter.Update()
    vtk_mesh = grad_filter.GetOutput()
    
    # Compute gradient of 'ab' (apicobasal)
    grad_filter2 = vtk.vtkGradientFilter()
    grad_filter2.SetInputData(vtk_mesh)
    grad_filter2.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "ab")
    grad_filter2.SetResultArrayName("grad_ab")
    grad_filter2.Update()
    vtk_mesh = grad_filter2.GetOutput()
    
    # Get arrays as numpy
    tm = numpy_support.vtk_to_numpy(vtk_mesh.GetPointData().GetArray("tm"))
    grad_tm = numpy_support.vtk_to_numpy(vtk_mesh.GetPointData().GetArray("grad_tm"))
    grad_ab = numpy_support.vtk_to_numpy(vtk_mesh.GetPointData().GetArray("grad_ab"))
    points = numpy_support.vtk_to_numpy(vtk_mesh.GetPoints().GetData())
    
    num_points = points.shape[0]
    
    # We will compute the metric tensor for each point
    metric_data = np.zeros((num_points, 9))
    fiber_data = np.zeros((num_points, 3))
    
    alpha_endo_rad = np.deg2rad(alpha_endo)
    alpha_epi_rad = np.deg2rad(alpha_epi)
    
    # Eigendecomposition of the metric tensor:
    # M = R * Lambda * R^T
    Lambda = np.diag([1.0 / (h_f**2), 1.0 / (h_t**2), 1.0 / (h_t**2)])
    
    print("Computing local basis, fibers, and metric tensors...")
    for i in range(num_points):
        g_tm = grad_tm[i]
        g_ab = grad_ab[i]
        
        # Normalize e_tm
        norm_tm = np.linalg.norm(g_tm)
        if norm_tm > 1e-12:
            e_tm = g_tm / norm_tm
        else:
            e_tm = np.array([1.0, 0.0, 0.0])
            
        # Orthogonalize e_ab against e_tm
        e_ab = g_ab - np.dot(g_ab, e_tm) * e_tm
        norm_ab = np.linalg.norm(e_ab)
        if norm_ab > 1e-12:
            e_ab = e_ab / norm_ab
        else:
            tmp = np.array([0.0, 1.0, 0.0])
            e_ab = tmp - np.dot(tmp, e_tm) * e_tm
            e_ab = e_ab / np.linalg.norm(e_ab)
            
        # e_rt = e_tm x e_ab
        e_rt = np.cross(e_tm, e_ab)
        
        # Rule-based angles
        t = tm[i]
        t = max(0.0, min(1.0, t))
        alpha = alpha_endo_rad * (1.0 - t) + alpha_epi_rad * t
        
        # Fiber, sheet, normal
        f = e_rt * np.cos(alpha) + e_ab * np.sin(alpha)
        f = f / np.linalg.norm(f) # Re-normalize to be safe
        s = e_tm
        n = np.cross(f, s)
        n = n / np.linalg.norm(n)
        
        R = np.column_stack((f, s, n))
        M = R @ Lambda @ R.T
        
        metric_data[i] = M.flatten()
        fiber_data[i] = f
        
    print(f"Writing metric tensor field to {output_msh}...")
    in_mesh = meshio.read(input_vtu)
    
    with open(output_msh, "w") as f:
        f.write('View "metric" {\n')
        for cell_block in in_mesh.cells:
            if cell_block.type == "tetra":
                cells = cell_block.data
                for cell in cells:
                    pts = points[cell] # 4x3
                    ms = metric_data[cell] # 4x9
                    p_str = ",".join([f"{p[0]},{p[1]},{p[2]}" for p in pts])
                    m_str = ",".join([",".join(map(str, m)) for m in ms])
                    f.write(f"ST({p_str}){{{m_str}}};\n")
        f.write("};\n")
        
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="meshes/idealizedBiv_cobiveco_result.vtu")
    parser.add_argument("--output", default="metric.pos")
    parser.add_argument("--hf", type=float, default=0.6, help="Longitudinal target edge length (mm)")
    parser.add_argument("--ht", type=float, default=0.15, help="Transverse target edge length (mm)")
    parser.add_argument("--alpha_endo", type=float, default=60.0, help="Helix angle at endocardium (deg)")
    parser.add_argument("--alpha_epi", type=float, default=-60.0, help="Helix angle at epicardium (deg)")
    args = parser.parse_args()
    
    generate_metric(args.input, args.output, args.hf, args.ht, args.alpha_endo, args.alpha_epi)
