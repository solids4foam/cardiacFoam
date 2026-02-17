import numpy as np
import pyvista as pv


"""
purkinje.py

Function to identify and tag elements with a Purkinje layer in the endocardial layer based on transmural distance.

Elements are tagged with int 1 if any of their nodes are within 10% transmural depth and belong to either the LV or RV.
"""



def add_purkinje_layer(
    mesh: pv.DataSet,
    transmural_min: float = 0.0,
    transmural_max: float = 0.1,
    lv_value: int = -1,
    rv_value: int = 1,
    field_name: str = "purkinjeLayer",
    purkinje_mult: float | None = None,
    diffusivity_field: str = "Diffusivity",
) -> pv.DataSet:
    """

    Parameters
    ----------
    mesh : pyvista.DataSet
        The input mesh with point data fields 'uvc_transmural' and 'uvc_intraventricular'.

    Returns
    -------
    mesh : pyvista.DataSet
        The same mesh, with an added cell data array 'purkinjeLayer' tagging
        endocardial volumes in the LV and RV. If a Diffusivity tensor exists
        and purkinje_mult is provided, the tensor is scaled in tagged cells.
    """
    if "uvc_transmural" not in mesh.point_data or "uvc_intraventricular" not in mesh.point_data:
        raise RuntimeError(
            "Missing required point_data fields 'uvc_transmural' and 'uvc_intraventricular'."
        )

    uvc_transmural = mesh.point_data["uvc_transmural"]
    uvc_intraventricular = mesh.point_data["uvc_intraventricular"]

    # Boolean masks
    is_endo = (uvc_transmural >= transmural_min) & (uvc_transmural <= transmural_max)
    is_rv = uvc_intraventricular == rv_value
    is_lv = uvc_intraventricular == lv_value

    purkinje_labels = np.where(is_endo & (is_rv | is_lv), 1, 0)

    # Cell-level tag
    cell_tags = np.zeros(mesh.n_cells, dtype=float)
    for i in range(mesh.n_cells):
        point_ids = mesh.get_cell(i).point_ids
        if np.any(purkinje_labels[point_ids] == 1):
            cell_tags[i] = 1

    mesh.cell_data[field_name] = cell_tags
    percent = (transmural_max - transmural_min) * 100.0
    print(
        f"Added purkinje layer in {percent:g}% of transmural depth "
        f"({np.count_nonzero(cell_tags)} cells tagged as {field_name})."
    )

    if purkinje_mult is not None and diffusivity_field in mesh.cell_data:
        diffusivity = mesh.cell_data[diffusivity_field]
        if diffusivity.ndim != 2 or diffusivity.shape[1] != 9:
            raise RuntimeError(
                f"Expected {diffusivity_field} to be (n_cells, 9), got {diffusivity.shape}."
            )
        mask = cell_tags == 1
        diffusivity = diffusivity.copy()
        diffusivity[mask] *= purkinje_mult
        mesh.cell_data[diffusivity_field] = diffusivity
        print(
            f"Scaled {diffusivity_field} by {purkinje_mult}x on "
            f"{np.count_nonzero(mask)} purkinjeLayer cells."
        )
    elif purkinje_mult is not None:
        print(
            f"Warning: {diffusivity_field} not found; purkinjeLayer tagged without scaling."
        )

    return mesh
