import numpy as np
import pyvista as pv


"""
purkinje.py

Function to identify and tag elements with a Purkinje layer in the endocardial layer based on transmural distance.

Elements are tagged with int 1 if any of their nodes are within 10% transmural depth and belong to either the LV or RV.
"""



def add_purkinje_layer(mesh: pv.DataSet) -> pv.DataSet:
    """

    Parameters
    ----------
    mesh : pyvista.DataSet
        The input mesh with point data fields 'uvc_transmural' and 'uvc_intraventricular'.

    Returns
    -------
    mesh : pyvista.DataSet
        The same mesh, with an added cell data array 'purkinjeLayer' tagging
        endocardial volumes in the LV and RV.
    """
    uvc_transmural = mesh.point_data["uvc_transmural"]
    uvc_intraventricular = mesh.point_data["uvc_intraventricular"]

    # Boolean masks
    is_endo = (uvc_transmural >= 0.0) & (uvc_transmural <= 0.1)
    is_rv = uvc_intraventricular == 1
    is_lv = uvc_intraventricular == -1

    purkinje_labels = np.where(is_endo & (is_rv | is_lv), 1, 0)

    # Cell-level tag
    cell_tags = np.zeros(mesh.n_cells, dtype=float)
    for i in range(mesh.n_cells):
        point_ids = mesh.get_cell(i).point_ids
        if np.any(purkinje_labels[point_ids] == 1):
            cell_tags[i] = 1

    mesh.cell_data["purkinjeLayer"] = cell_tags
    print(f"Number of mesh volumes tagged as purkinjeLayer = {np.count_nonzero(cell_tags)}")


    return mesh
    return mesh
