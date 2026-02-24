"""
scar.py

Helpers to tag scar cells and optionally zero diffusivity tensors.
"""


def add_scar_from_selection(
    full_mesh,
    selection_mesh,
    id_array: str = "GlobalCellIds",
    scar_name: str = "Scar",
    diffusivity_name: str = "Diffusivity",
    diffusivity_scale: float = 0.0,
    scar_value: float = 1.0,
    return_mask: bool = False,
):
    """
    Tag scar cells based on matching cell IDs from a selection mesh.

    Parameters
    ----------
    full_mesh : pv.DataSet
        Full mesh containing all cells.
    selection_mesh : pv.DataSet
        Selection mesh containing a subset of cells.
    id_array : str
        Cell data array name containing global cell IDs.
    scar_name : str
        Cell data array name to store the scar tag.
    diffusivity_name : str
        Cell data array name for the diffusivity tensor.
    diffusivity_scale : float
        Scale diffusivity tensor in scar cells (0.0 = zero, 1.0 = unchanged).
    scar_value : float
        Value assigned to scar cells (typically 0.0 to 1.0).
    return_mask : bool
        Return a boolean scar mask along with the mesh.
    """
    if id_array not in full_mesh.cell_data:
        raise RuntimeError(
            f"full_mesh has no CellData '{id_array}'. "
            "In ParaView, run 'Generate Global IDs' on the full mesh and save it."
        )

    if id_array not in selection_mesh.cell_data:
        raise RuntimeError(
            f"selection_mesh has no CellData '{id_array}'. "
            "In ParaView, run 'Generate Global IDs' BEFORE extracting selection."
        )

    import numpy as np

    full_ids = np.asarray(full_mesh.cell_data[id_array])
    sel_ids = np.unique(np.asarray(selection_mesh.cell_data[id_array]))

    scar_mask = np.isin(full_ids, sel_ids)
    scar = np.ones(full_ids.shape, dtype=np.float32)
    scar[scar_mask] = scar_value

    full_mesh.cell_data[scar_name] = scar

    if diffusivity_name in full_mesh.cell_data:
        diffusivity = np.array(full_mesh.cell_data[diffusivity_name], copy=True)
        diffusivity[scar_mask] = diffusivity[scar_mask] * diffusivity_scale
        full_mesh.cell_data[diffusivity_name] = diffusivity
    else:
        print(
            f"Warning: '{diffusivity_name}' not found in full mesh cell data; "
            "no tensor update applied."
        )

    if return_mask:
        return full_mesh, scar_mask
    return full_mesh
