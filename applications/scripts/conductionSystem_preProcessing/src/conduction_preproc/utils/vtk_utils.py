"""
vtk_utils.py

Utility functions for handling VTK mesh files and ASCII cleanup.

Functions:
    read_values(lines, i, nvals)
        Read a fixed number of values from lines of a VTK ASCII file.

    inspect_fields(infile)
        Print available point and cell data fields from a VTK file.

    interactive_remove_blank_lines(filepath, output_file=None)
        Interactively remove extra blank lines from an ASCII VTK file to make it valid for OpenFOAM reader.
"""

def read_values(lines, i, nvals):
    vals = []
    while i < len(lines) and len(vals) < nvals:
        vals.extend(lines[i].split())
        i += 1
    if len(vals) != nvals:
        raise ValueError(f"Expected {nvals} values, got {len(vals)}")
    return vals, i



def inspect_fields(infile):
    """Check FIELD, SCALARS, VECTORS, TENSORS inside a VTK file for debugging."""
    with open(infile, "r") as f:
        lines = f.readlines()

    current_section = None
    for idx, line in enumerate(lines):
        stripped = line.strip()

        # Track current section
        if stripped.startswith("POINT_DATA") or stripped.startswith("CELL_DATA"):
            current_section = stripped.split()[0]
            print(f"\nSection {current_section} at line {idx}")
            continue

        # FIELD arrays
        if stripped.startswith("FIELD") and current_section:
            tokens = stripped.split()
            if len(tokens) < 3:
                continue
            field_name = tokens[1]
            n_arrays = int(tokens[2])
            print(f"  FIELD '{field_name}' starting at line {idx} with {n_arrays} arrays")

            i = idx + 1
            for arr_idx in range(n_arrays):
                if i >= len(lines):
                    break
                hdr = lines[i].strip()
                parts = hdr.split()
                if len(parts) != 4:
                    print(f"    Malformed array header at line {i}: {hdr}")
                    i += 1
                    continue

                name, ncomp_s, ntup_s, dtype = parts
                try:
                    ncomp = int(ncomp_s)
                    ntup = int(ntup_s)
                except ValueError:
                    print(f"    Invalid nComp or nTuples at line {i}: {hdr}")
                    i += 1
                    continue

                print(f"    Array #{arr_idx+1}: {name}, nComp={ncomp}, nTuples={ntup}, type={dtype}")
                i += 1
                nvals = ncomp * ntup
                _, i = read_values(lines, i, nvals)

        # Scalars, Vectors, Tensors
        if current_section and (stripped.startswith("SCALARS") or
                                stripped.startswith("VECTORS") or
                                stripped.startswith("TENSORS")):
            tokens = stripped.split()
            arr_type = tokens[0]
            name = tokens[1] if len(tokens) > 1 else "<unnamed>"
            dtype = tokens[2] if len(tokens) > 2 else "<unknown>"
            print(f"  {arr_type} '{name}' ({dtype}) already in place at line {idx}")



def interactive_remove_blank_lines(filepath, output_file=None):
    """Remove selected blank lines from a file (interactive)."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    blank_map = {}
    count = 0
    for i, line in enumerate(lines):
        if line.strip() == "":
            before = lines[i - 1].rstrip("\n") if i > 0 else "<START>"
            after = lines[i + 1].rstrip("\n") if i < len(lines) - 1 else "<END>"

            print(f"\n--- Blank line {count} (at file line {i + 1}) ---")
            print(f"Before: {before}")
            print(f"After : {after}")

            blank_map[count] = i
            count += 1

    if not blank_map:
        print("No blank lines found.")
        return

    choice = input("\nEnter blank line counts to remove (comma-separated), or 'all': ").strip()
    if choice.lower() == "all":
        indices_to_remove = list(blank_map.values())
    elif choice:
        ids = [int(x) for x in choice.split(",") if x.strip().isdigit()]
        indices_to_remove = [blank_map[i] for i in ids if i in blank_map]
    else:
        print("No lines selected. Nothing removed.")
        return

    new_lines = [line for i, line in enumerate(lines) if i not in indices_to_remove]
    if output_file is None:
        output_file = filepath
    with open(output_file, "w") as f:
        f.writelines(new_lines)

    print(f"\nRemoved {len(indices_to_remove)} blank lines. Output written to '{output_file}'.")

def remove_blank_lines(filepath, output_file=None):
    """Remove all blank lines from a file (non-interactive)."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    new_lines = [line for line in lines if line.strip() != ""]
    if output_file is None:
        output_file = filepath
    with open(output_file, "w") as f:
        f.writelines(new_lines)

    print(f"\nRemoved {len(lines) - len(new_lines)} blank lines. Output written to '{output_file}'.")


def plot_mesh_surface(mesh, scalar_name: str = None, opacity: float = 0.5, cmap: str = "viridis", save_screenshot: str = None):
    """
    Plot only the surface of a mesh for faster visualization.

    Parameters
    ----------
    mesh : pyvista.DataSet
        The mesh to visualize.
    scalar_name : str, optional
        Name of a cell or point data array to color by.
    opacity : float, optional
        Mesh opacity (0.0 - 1.0).
    cmap : str, optional
        Colormap name.
    save_screenshot : str, optional
        File path to save a screenshot instead of showing interactively.
    """
    import pyvista as pv

    # Extract only the outer surface
    surface = mesh.extract_surface()

    plotter = pv.Plotter()
    plotter.add_mesh(
        surface,
        scalars=scalar_name,
        cmap=cmap,
        opacity=opacity,
        show_edges=True
    )

    if save_screenshot:
        plotter.screenshot(save_screenshot)
        plotter.close()
        print(f"Screenshot saved to {save_screenshot}")
    else:
        plotter.show()

