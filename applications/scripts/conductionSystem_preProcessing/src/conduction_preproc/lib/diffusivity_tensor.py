"""
tensors.py

Functions to compute and attach diffusivity tensors to a cardiac mesh.

Functions:
    normalize(v)
        Normalize a 3-component vector. Returns a zero vector if norm is zero.

    add_diffusivity_tensor_ventricles(mesh, D_f, D_s, D_n,
                                       fiber_name="fiber", sheet_name="sheet",
                                       purkinje_name="purkinjeLayer",
                                       tensor_name="Diffusivity")
        Add a 3x3 diffusion tensor for each cell in the mesh based on the local
        fiber, sheet, and Purkinje layer information. 
        Can be extended for atria.
"""

import numpy as np


def normalize(v):
    """
    Normalize a 3D vector.
    """
    norm = np.linalg.norm(v)
    return v / norm if norm > 0 else np.zeros_like(v)


def add_diffusivity_tensor_ventricles(mesh, D_f, D_s, D_n,
                                      fiber_name="fiber", sheet_name="sheet",
                                      purkinje_name="purkinjeLayer",
                                      tensor_name="Diffusivity",
                                      purkinje_mult=2.0):
    """
    Compute and attach diffusion tensors to each cell in a mesh.
    Parameters
    ----------
    mesh : pyvista.DataSet
        The mesh containing cell data fields for fibers, sheets, and Purkinje tags.
    D_f : float
        Diffusivity along fiber direction.
    D_s : float
        Diffusivity along sheet direction.
    D_n : float
        Diffusivity along normal direction.
    fiber_name : str, optional
        Name of the fiber vector field in cell_data.
    sheet_name : str, optional
        Name of the sheet vector field in cell_data.
    purkinje_name : str, optional
        Name of the Purkinje layer scalar field in cell_data.
    tensor_name : str, optional
        Name to assign the new diffusion tensor field.

    Returns
    -------
    mesh : pyvista.DataSet
        Mesh with new cell_data[tensor_name] added.
    """
    print("\n-------------------------\nAdding Diffusivity Tensor\n-------------------------")

    fiber_vecs = mesh.cell_data[fiber_name]
    sheet_vecs = mesh.cell_data[sheet_name]
    if purkinje_name in mesh.cell_data:
        purkinje_scalars = mesh.cell_data[purkinje_name]
    else:
        purkinje_scalars = np.zeros(mesh.n_cells, dtype=fiber_vecs.dtype)

    # Normalize vectors
    fiber_vecs = np.array([normalize(f) for f in fiber_vecs])
    sheet_vecs = np.array([normalize(s) for s in sheet_vecs])
    normal_vecs = np.array([normalize(np.cross(f, s)) for f, s in zip(fiber_vecs, sheet_vecs)])

    # Compute tensors
    D_tensors = []
    for f, s, n, p in zip(fiber_vecs, sheet_vecs, normal_vecs, purkinje_scalars):
        scale = 1 + (purkinje_mult - 1) * p
        D = scale * (
            D_f * np.outer(f, f) + D_s * np.outer(s, s) + D_n * np.outer(n, n)
        )
        D_tensors.append(D.flatten())

    D_tensors = np.array(D_tensors, dtype=fiber_vecs.dtype)
    mesh.cell_data[tensor_name] = D_tensors
    print(f"Added tensor '{tensor_name}' shape: {D_tensors.shape}")
    return mesh
