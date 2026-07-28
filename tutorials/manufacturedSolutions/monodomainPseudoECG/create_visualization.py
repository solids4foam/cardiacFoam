import pyvista as pv
import os

def create_vtk():
    foam_file = "monodomainPseudoECG.foam"
    if not os.path.exists(foam_file):
        with open(foam_file, 'w') as f:
            pass

    # Read OpenFOAM file
    print("Reading OpenFOAM data...")
    reader = pv.OpenFOAMReader(foam_file)
    time_values = reader.time_values
    if time_values:
        print(f"Setting time to {time_values[-1]}")
        reader.set_active_time_value(time_values[-1])

    mesh = reader.read()

    # Extract internal mesh
    if "internalMesh" in mesh.keys():
        internal_mesh = mesh["internalMesh"]
    else:
        # Fallback if no internalMesh name
        internal_mesh = mesh[0]

    # Probe coordinates
    probes = {
        "E1": (-0.5, 0.5, 0.5),
        "E2": (1.5, 0.5, 0.5),
        "E3": (1.2, 0.23, 0.61),
        "E4": (1.35, 0.74, 0.28),
        "E5": (1.55, 0.41, 0.83)
    }

    # We will combine everything into a MultiBlock dataset
    output = pv.MultiBlock()
    output.append(internal_mesh, "Cube_Voltage")

    # Add spheres for probes
    for name, pos in probes.items():
        sphere = pv.Sphere(radius=0.05, center=pos)
        output.append(sphere, f"Probe_{name}")

    output_filename = "voltage_and_probes.vtm"
    output.save(output_filename)
    print(f"Successfully saved {output_filename}")

if __name__ == "__main__":
    create_vtk()
