"""
Prepare WebGL-ready simulation data from VTP boundary patches.

Reads all VTP surface patches from the realBiv02Case, decimates each individually,
concatenates into a single non-indexed geometry with per-patch explosion vectors
and scalar color buffers.

Usage:
    python prepare_simulation_data.py <patches_dir>
    
Example:
    python prepare_simulation_data.py /path/to/realBiv02Case/patches
"""
import sys
import json
import os
import warnings

warnings.filterwarnings("ignore")

try:
    import pyvista as pv
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError:
    print("PyVista, NumPy, and Matplotlib are required.")
    sys.exit(1)

# Per-patch decimation targets: large patches get more aggressive reduction
# to keep total vertex count manageable for WebGL
DECIMATION = {
    "epi":      0.95,    # 142k faces → ~7k
    "endo_lv":  0.95,    # 76k faces → ~4k
    "endo_rv":  0.95,    # 80k faces → ~4k
    "tag_3":    0.95,    # 70k faces → ~3.5k
    "tag_4":    0.95,    # 50k faces → ~2.5k
    "tag_5":    0.92,    # 21k faces → ~1.7k
    "tag_6":    0.90,    # 7k faces → ~750
    # Smaller patches: lighter decimation or none
}
DEFAULT_DECIMATION = 0.80  # For small patches (<5k faces)
MIN_FACES_TO_DECIMATE = 500  # Don't decimate patches smaller than this


def main():
    if len(sys.argv) < 2:
        print("Usage: python prepare_simulation_data.py <patches_dir>")
        sys.exit(1)
    
    patches_dir = sys.argv[1]
    out_dir = "frontend/public/scalars"
    os.makedirs(out_dir, exist_ok=True)
    
    # Load the manifest
    manifest_path = os.path.join(patches_dir, "all_patches_manifest.json")
    if not os.path.exists(manifest_path):
        print(f"ERROR: {manifest_path} not found. Run extract_all_surfaces.py first.")
        sys.exit(1)
    
    with open(manifest_path) as f:
        patch_manifest = json.load(f)
    
    print(f"Loading {len(patch_manifest)} patches from {patches_dir} ...")
    
    # ── Load + Decimate each patch ──
    all_positions = []    # List of (N, 3) float32 arrays
    all_patch_ids = []    # Patch ID per expanded vertex (for explosion)
    all_scalars = {}      # field_name → list of arrays
    
    patch_centers = {}    # patch_name → center point (for explosion vectors)
    total_expanded = 0
    
    # Collect all field names across patches
    field_names = set()
    
    for patch_name, info in sorted(patch_manifest.items()):
        vtp_path = os.path.join(patches_dir, info["file"])
        if not os.path.exists(vtp_path):
            print(f"  SKIP {patch_name}: {vtp_path} not found")
            continue
        
        surf = pv.read(vtp_path)
        if not surf.is_all_triangles:
            surf = surf.triangulate()
        
        orig_faces = surf.n_cells
        
        # Convert cell data to point data BEFORE decimation
        if surf.cell_data:
            surf = surf.cell_data_to_point_data()
        
        # ── Save scalar fields BEFORE decimation (decimate strips point_data!) ──
        pre_dec_points = surf.points.copy()
        pre_dec_scalars = {}
        for k in list(surf.point_data.keys()):
            if k.startswith("vtk") or k in ("RGB", "Normals", "patchCode", "patchId",
                                              "patchName", "sourceId", "fiber", "sheet"):
                continue
            pre_dec_scalars[k] = np.asarray(surf.point_data[k], dtype=np.float32)
            field_names.add(k)
        
        # Decimate
        dec_ratio = DECIMATION.get(patch_name, DEFAULT_DECIMATION)
        if orig_faces > MIN_FACES_TO_DECIMATE:
            surf = surf.decimate(dec_ratio)
        
        n_faces = surf.n_cells
        
        # KDTree remap: map decimated vertices back to pre-decimated scalar values
        from scipy.spatial import cKDTree
        tree = cKDTree(pre_dec_points)
        _, nn_idx = tree.query(surf.points)
        
        # Expand to non-indexed (per-face vertices for clean explosion)
        faces = surf.faces.reshape(-1, 4)[:, 1:].astype(int)
        expanded_pts = surf.points[faces.flatten()].reshape(-1, 3).astype(np.float32)
        expanded_nn_idx = nn_idx[faces.flatten()]
        n_expanded = len(expanded_pts)
        
        # Record center for explosion vectors
        patch_centers[patch_name] = np.mean(expanded_pts, axis=0)
        
        # Store positions
        all_positions.append(expanded_pts)
        
        # Store patch ID for every expanded vertex (for explosion grouping)
        all_patch_ids.append(np.full(n_expanded, len(all_positions) - 1, dtype=np.int32))
        
        # Collect remapped scalar fields
        for k, pre_arr in pre_dec_scalars.items():
            expanded_arr = pre_arr[expanded_nn_idx]
            if k not in all_scalars:
                all_scalars[k] = []
            all_scalars[k].append(expanded_arr)
        
        total_expanded += n_expanded
        print(f"  {patch_name:15s}: {orig_faces:6d} → {n_faces:5d} faces,"
              f" {n_expanded:6d} expanded verts ({dec_ratio*100:.0f}% dec)"
              f"  [{len(pre_dec_scalars)} fields]")
    
    # ── Concatenate all patches ──
    print(f"\nConcatenating {len(all_positions)} patches → {total_expanded} total expanded verts ...")
    
    combined_pts = np.vstack(all_positions)
    combined_patch_ids = np.concatenate(all_patch_ids)
    
    # Center at origin
    center = combined_pts.mean(axis=0)
    combined_pts -= center
    
    # ── Rotate so apex-base axis aligns with Y-up ──
    # Apex-base direction computed from UVC longitudinal correlation with XYZ
    # (base → apex direction: negative X, positive Y, slight Z)
    apex_base = np.array([-0.63192916, 0.7414862, 0.22552988])
    apex_base /= np.linalg.norm(apex_base)
    target_up = np.array([0.0, 1.0, 0.0])
    
    # Rotation matrix via Rodrigues formula: rotate apex_base → target_up
    v = np.cross(apex_base, target_up)
    s = np.linalg.norm(v)
    c = np.dot(apex_base, target_up)
    if s > 1e-6:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * (1 - c) / (s * s)
        combined_pts = (R @ combined_pts.T).T.astype(np.float32)
        print(f"  Rotated apex-base to Y-up (rotation angle: {np.degrees(np.arccos(c)):.1f}°)")
    
    # Recompute patch centers relative to new origin + rotation
    for k in patch_centers:
        pc = patch_centers[k] - center
        if s > 1e-6:
            pc = R @ pc
        patch_centers[k] = pc
    
    # ── Two-Level Hierarchical Explosion ──
    # Phase 1 (macro): chambers separate from global center (all sub-parts move together)
    # Phase 2 (micro): sub-parts within each chamber separate from chamber barycenter
    print("Computing hierarchical explosion vectors ...")
    global_center = np.mean(combined_pts, axis=0)
    
    # Load chamber groupings from patch_colors.json
    color_map_path = os.path.join(out_dir, "patch_colors.json")
    chamber_groups = {}  # chamber_name -> list of patch names
    if os.path.exists(color_map_path):
        with open(color_map_path) as f:
            cm = json.load(f)
        for chamber_name, ch_data in cm.get("chambers", {}).items():
            chamber_groups[chamber_name] = ch_data["patches"]
    
    # Build patch_name -> chamber_name lookup
    patch_to_chamber = {}
    for chamber_name, patch_list in chamber_groups.items():
        for pname in patch_list:
            patch_to_chamber[pname] = chamber_name
    
    # Build sorted patch list with indices
    sorted_patches = sorted(patch_manifest.items())
    
    # Compute chamber barycenters (average of all vertices belonging to all patches in chamber)
    chamber_centers = {}
    for chamber_name in chamber_groups:
        chamber_mask = np.zeros(len(combined_pts), dtype=bool)
        for i, (pname, _) in enumerate(sorted_patches):
            if patch_to_chamber.get(pname) == chamber_name:
                chamber_mask |= (combined_patch_ids == i)
        if np.any(chamber_mask):
            chamber_centers[chamber_name] = np.mean(combined_pts[chamber_mask], axis=0)
    
    # Phase 1: macro vectors — each vertex moves with its chamber away from global center
    macro_vecs = np.zeros_like(combined_pts)
    for chamber_name, ch_center in chamber_centers.items():
        for i, (pname, _) in enumerate(sorted_patches):
            if patch_to_chamber.get(pname) == chamber_name:
                mask = (combined_patch_ids == i)
                macro_vecs[mask] = ch_center - global_center
    
    # Phase 2: micro vectors — each patch moves away from its chamber center
    micro_vecs = np.zeros_like(combined_pts)
    for i, (pname, _) in enumerate(sorted_patches):
        chamber_name = patch_to_chamber.get(pname)
        if chamber_name and chamber_name in chamber_centers:
            mask = (combined_patch_ids == i)
            if np.any(mask):
                patch_center = np.mean(combined_pts[mask], axis=0)
                micro_vecs[mask] = patch_center - chamber_centers[chamber_name]
    
    out_macro = os.path.join(out_dir, "disassemble_macro.bin")
    macro_vecs.astype(np.float32).flatten().tofile(out_macro)
    print(f"  → disassemble_macro.bin ({os.path.getsize(out_macro)/1024:.1f} KB)")
    
    out_micro = os.path.join(out_dir, "disassemble_micro.bin")
    micro_vecs.astype(np.float32).flatten().tofile(out_micro)
    print(f"  → disassemble_micro.bin ({os.path.getsize(out_micro)/1024:.1f} KB)")
    
    # ── JSON Geometry ──
    geom_data = {
        "positions": [round(float(x), 3) for x in combined_pts.flatten()],
        "vertexCount": int(len(combined_pts)),
        "faceCount": int(len(combined_pts) // 3)
    }
    geom_path = "frontend/public/heart_geometry.json"
    with open(geom_path, 'w') as f:
        json.dump(geom_data, f)
    print(f"  → heart_geometry.json ({os.path.getsize(geom_path)/1024:.0f} KB)")
    
    # ── Color Buffers ──
    web_manifest = {
        "vertices": int(len(combined_pts)),
        "faces": int(len(combined_pts) // 3),
        "fields": {}
    }
    cmap_viridis = plt.get_cmap("viridis")
    
    # Load curated anatomical color palette
    color_map_path = os.path.join(out_dir, "patch_colors.json")
    patch_colors = {}
    if os.path.exists(color_map_path):
        with open(color_map_path) as f:
            cm = json.load(f)
        patch_colors = {k: v["color"] for k, v in cm.get("patches", {}).items()}
        print(f"  Loaded {len(patch_colors)} curated colors from patch_colors.json")
    
    # Generate tags color buffer from curated per-patch colors
    tag_rgb = np.full((len(combined_pts), 3), 180, dtype=np.uint8)  # default light gray
    for i, (patch_name, _) in enumerate(sorted(patch_manifest.items())):
        mask = (combined_patch_ids == i)
        if np.any(mask) and patch_name in patch_colors:
            tag_rgb[mask] = patch_colors[patch_name]
    
    out_bin = os.path.join(out_dir, "tags.bin")
    tag_rgb.flatten().tofile(out_bin)
    n_patches = len(patch_manifest)
    web_manifest["fields"]["tags"] = {"min": 0, "max": n_patches, "file": "/scalars/tags.bin"}
    print(f"  → tags.bin ({os.path.getsize(out_bin)/1024:.1f} KB)")
    
    # Generate color buffers for scalar fields
    for field_name in sorted(field_names):
        if field_name not in all_scalars:
            continue
        
        scalars = np.concatenate(all_scalars[field_name])
        
        # For vector fields, take magnitude
        if len(scalars.shape) > 1 and scalars.shape[1] > 1:
            scalars = np.linalg.norm(scalars, axis=1)
        
        # Outlier culling (UVC atria defaults at 100.0 = null/undefined)
        valid = scalars < 99.0
        
        # UVC coordinates are defined in [0, 1] — use fixed range, not auto-scaled
        if field_name.startswith("uvc_"):
            lo, hi = 0.0, 1.0
        elif np.any(valid):
            lo, hi = float(np.nanmin(scalars[valid])), float(np.nanmax(scalars[valid]))
        else:
            lo, hi = 0.0, 1.0
        
        # Clamp: values ≥ 99 are null → set to 0 before normalization
        scalars = np.where(valid, scalars, 0.0)
        norm = (scalars - lo) / (hi - lo) if hi > lo else np.zeros_like(scalars)
        norm = np.clip(norm, 0, 1)
        rgb = (cmap_viridis(norm)[:, :3] * 255).astype(np.uint8)
        
        out_bin = os.path.join(out_dir, f"{field_name}.bin")
        rgb.flatten().tofile(out_bin)
        web_manifest["fields"][field_name] = {"min": lo, "max": hi, "file": f"/scalars/{field_name}.bin"}
        print(f"  → {field_name}.bin ({os.path.getsize(out_bin)/1024:.1f} KB)")
    
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(web_manifest, f, indent=2)
    
    total_kb = sum(os.path.getsize(os.path.join(out_dir, fn)) for fn in os.listdir(out_dir)) / 1024
    ply_kb = os.path.getsize(geom_path) / 1024
    print(f"\n✓ Payload: {total_kb:.0f} KB (scalars) + {ply_kb:.0f} KB (geometry)")
    print(f"✓ {len(combined_pts)} expanded verts | {len(combined_pts)//3} faces | {len(patch_manifest)} patches")


if __name__ == "__main__":
    main()
