import os
import json
import pyvista as pv
from multiprocessing import Pool, cpu_count

def process_timestep(args):
    i, t, case_path, purkinje_dir, closest_purkinje, out_dir = args
    
    # Setup OpenFOAM reader just for this time
    reader = pv.OpenFOAMReader(case_path)
    reader.enable_all_cell_arrays()
    reader.enable_all_point_arrays()
    reader.set_active_time_value(t)
    foam_data = reader.read()
    
    heart_mesh = foam_data["internalMesh"] if "internalMesh" in foam_data else foam_data[0]
    
    # 1. Extract surface and decimate!
    # VTK doesn't have GPU decimation, but we can do it in parallel on the CPU
    surface = heart_mesh.extract_surface(algorithm='dataset_surface')
    n_faces = getattr(surface, "n_cells", 0)
    target_faces = 150000
    
    if n_faces > target_faces:
        ratio = 1.0 - (target_faces / n_faces)
        # Using decimate_pro which is VTK's fast CPU decimator
        surface = surface.decimate_pro(
            ratio, feature_angle=15.0, splitting=True, 
            boundary_vertex_deletion=False, preserve_topology=False
        )
        
    # 2. Load Purkinje
    purkinje_vtk_path = os.path.join(purkinje_dir, closest_purkinje["name"])
    try:
        purkinje_mesh = pv.read(purkinje_vtk_path)
    except Exception as e:
        purkinje_mesh = pv.PolyData()
        
    import numpy as np
    
    # 3. Rename fields to match perfectly so PyVista merge preserves them!
    if "Vm_V" in purkinje_mesh.point_data:
        purkinje_mesh.point_data["Vm"] = purkinje_mesh.point_data.pop("Vm_V")
        
    # Add a RegionId field just in case we need to separate them visually later
    surface.point_data["RegionId"] = np.zeros(surface.n_points, dtype=np.int32)
    if purkinje_mesh.n_points > 0:
        purkinje_mesh.point_data["RegionId"] = np.ones(purkinje_mesh.n_points, dtype=np.int32)

    # 4. Save as a single merged PolyData
    combined = surface.merge(purkinje_mesh)
    
    filename = f"combined_step_{i:03d}.vtp"
    filepath = os.path.join(out_dir, filename)
    combined.save(filepath)
    
    print(f"Done step {i} (t={t:.4f}) -> {filename} (Faces: {getattr(combined, 'n_cells', 0)})")
    return {"name": filename, "time": float(t)}

def export_combined_data():
    case_path = "monodomainHeart.foam"
    purkinje_series_path = "postProcessing/purkinjeNetworkVTK/purkinjeNetwork.vtk.series"
    out_dir = "exported_combined_purkinje"
    series_json_path = os.path.join(out_dir, "combined_purkinje.vtk.series")

    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(purkinje_series_path):
        print(f"Error: {purkinje_series_path} not found.")
        return
        
    with open(purkinje_series_path, "r") as f:
        purkinje_series = json.load(f)
    
    purkinje_files = purkinje_series.get("files", [])
    purkinje_dir = os.path.dirname(purkinje_series_path)

    # Get time values
    reader = pv.OpenFOAMReader(case_path)
    time_values = reader.time_values
    
    print(f"Found {len(time_values)} time steps in OpenFOAM data.")
    print(f"Found {len(purkinje_files)} time steps in Purkinje network data.")
    print(f"Using 6 CPU cores to process all {len(time_values)} timesteps in parallel...\n")

    # Build arguments for the multiprocessing pool
    pool_args = []
    for i, t in enumerate(time_values):
        closest_purkinje = min(purkinje_files, key=lambda x: abs(x["time"] - t))
        pool_args.append((i, t, case_path, purkinje_dir, closest_purkinje, out_dir))

    # Process all timesteps in parallel using 6 CPU cores
    series_files = []
    with Pool(processes=6) as pool:
        for result in pool.imap(process_timestep, pool_args):
            series_files.append(result)

    # Sort files by time just to be safe
    series_files.sort(key=lambda x: x["time"])

    # Write the .vtk.series index file
    with open(series_json_path, "w") as f:
        json.dump({"file-series-version": "1.0", "files": series_files}, f, indent=4)
        
    print(f"\nDone! Exported full, optimized simulation data to: {series_json_path}")

if __name__ == "__main__":
    export_combined_data()
