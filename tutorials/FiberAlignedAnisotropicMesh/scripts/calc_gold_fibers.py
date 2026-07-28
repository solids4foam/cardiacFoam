import vtk
import numpy as np
import sys
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

def get_gradient(mesh, scalar_name, grad_name):
    grad_filter = vtk.vtkGradientFilter()
    grad_filter.SetInputData(mesh)
    grad_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, scalar_name)
    grad_filter.SetResultArrayName(grad_name)
    grad_filter.Update()
    return grad_filter.GetOutput()

def normalize(v):
    norm = np.linalg.norm(v, axis=1, keepdims=True)
    norm[norm == 0] = 1e-12
    return v / norm

def main():
    in_vtu = sys.argv[1]
    out_vtu = sys.argv[2]

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(in_vtu)
    reader.Update()
    mesh = reader.GetOutput()

    # Compute grad(tm) -> eT raw
    mesh = get_gradient(mesh, 'tm', 'grad_tm')
    # Compute grad(ab) -> eL raw
    mesh = get_gradient(mesh, 'ab', 'grad_ab')

    pd = mesh.GetPointData()
    grad_tm = vtk_to_numpy(pd.GetArray('grad_tm'))
    grad_ab = vtk_to_numpy(pd.GetArray('grad_ab'))
    tm = vtk_to_numpy(pd.GetArray('tm'))

    # Orthonormalization (Bayer/Cobiveco convention)
    eT = normalize(grad_tm)
    
    # Gram-Schmidt for eL
    eL_raw = grad_ab
    dot_ab_eT = np.sum(eL_raw * eT, axis=1, keepdims=True)
    eL = normalize(eL_raw - dot_ab_eT * eT)
    
    # eC is cross(eL, eT)
    eC = np.cross(eL, eT)
    
    # Rotation angle alpha based on tm.
    # In Cobiveco, tm is 0 at Epi, 1 at Endo usually. Wait, in idealizedBiv_cobiveco_result.vtu, let's check tm bounds.
    # We will just linearly interpolate from -60 (min tm) to +60 (max tm)
    tm_min = np.min(tm)
    tm_max = np.max(tm)
    # Normalize tm to 0-1 range
    tm_norm = (tm - tm_min) / (tm_max - tm_min)
    
    alpha_deg = (-60.0 * (1.0 - tm_norm)) + (60.0 * tm_norm)
    alpha_rad = np.radians(alpha_deg)[:, np.newaxis]
    
    # Rotate eC around eT by alpha: f0 = cos(alpha)*eC + sin(alpha)*eL
    f0 = np.cos(alpha_rad) * eC + np.sin(alpha_rad) * eL
    f0 = normalize(f0)
    
    # Add f0 back to VTU
    f0_vtk = numpy_to_vtk(f0, deep=True)
    f0_vtk.SetName("f0")
    pd.AddArray(f0_vtk)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(out_vtu)
    writer.SetInputData(mesh)
    writer.Write()
    
    print(f"Calculated f0 and wrote to {out_vtu}")

if __name__ == "__main__":
    main()
