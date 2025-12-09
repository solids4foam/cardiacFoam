
"""
main.py

Orchestrates mesh processing:
using the config.py settings, it performs the following steps:
1. Inspect fields in input VTK.
2. Optionally Purkinje layer.
3. Optionally add diffusivity tensors.
4. Convert FIELD arrays to SCALARS/VECTORS/TENSORS for OpenFOAM.
"""
import os
import shutil
import pyvista as pv

from config import INPUT_FILE, MIDDLE_FILE, OUTPUT_FILE, Ventricular_DF, Ventricular_DS, Ventricular_DN
from lib.purkinje import add_purkinje_layer
from lib.diffusivity_tensor import add_diffusivity_tensor_ventricles
from utils.vtk_utils import inspect_fields, interactive_remove_blank_lines, plot_mesh_surface
from utils.vtk_convert_arrays_to_fields import convert_vtk_file






os.makedirs("output", exist_ok=True)



def main():
    # Step 1: Inspect input mesh
    print(f"Loading mesh: {INPUT_FILE}")
    inspect_fields(INPUT_FILE)

    # Step 2: Add Purkinje layer
    mesh = pv.read(INPUT_FILE)

    # --- Prompt for Purkinje layer ---
    add_purkinje = input("\nDo you want to add the Purkinje layer? (y/n): ").strip().lower()
    if add_purkinje == "y":
        mesh = add_purkinje_layer(mesh)
        mesh.save(MIDDLE_FILE, binary=False)
        print(f"Purkinje layer added. Mesh saved to {MIDDLE_FILE}")
        
    else:
        # If skipping, just copy input to middle output
        shutil.copyfile(INPUT_FILE, MIDDLE_FILE)
        print(f"Purkinje layer skipped. Mesh copied to {MIDDLE_FILE}")
    # Plot only the surface of the Purkinje layer
    



    inspect_fields(MIDDLE_FILE)

    # Step 4: Optionally add diffusivity tensor
    add_diff = input("\nDo you want to add a diffusivity tensor to the mesh? (y/n): ").strip().lower()
    if add_diff == "y":
        mesh = add_diffusivity_tensor_ventricles(mesh, Ventricular_DF, Ventricular_DS, Ventricular_DN)
        mesh.save(OUTPUT_FILE, binary=False)
        print(f"Diffusivity tensor added. Mesh saved to {OUTPUT_FILE}")

        # Clean blank lines in output
        interactive_remove_blank_lines(OUTPUT_FILE, OUTPUT_FILE)
    else:
        # Just copy middle mesh to final output
        shutil.copyfile(MIDDLE_FILE, OUTPUT_FILE)
        print(f"No tensor added. Mesh copied to {OUTPUT_FILE}")

    # Step 5: Convert FIELD arrays to SCALARS/VECTORS/TENSORS
    convert_vtk_file(OUTPUT_FILE, OUTPUT_FILE)
    print("VTK FIELD arrays converted for readable OpenFOAM format: SCALARS/VECTORS/TENSORS Fields.")
    inspect_fields(OUTPUT_FILE)
    print(f"\nDone 😀. Please verify the fields and type.")



if __name__ == "__main__":
    main()
