import numpy as np
import pyvista as pv
from pathlib import Path

def rigid_transform_3D(A, B):
    """ Rigid transformation mapping A -> B (Kabsch Algorithm) """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(B)):
        raise ValueError("Landmarks contain non-finite values (NaN/Inf). Cannot align.")
        
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    H = (A - centroid_A).T @ (B - centroid_B)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        print("  Handedness mismatch detected; adjusting rotation.")
        Vt[2, :] *= -1
        R = Vt.T @ U.T
    return R, centroid_B - R @ centroid_A

def rotate_around_axis(points, p1, p2, angle_deg):
    """ Rotates points around the axis defined by p1 and p2 by angle_deg """
    angle_rad = np.radians(angle_deg)
    axis = p2 - p1
    dist = np.linalg.norm(axis)
    if dist < 1e-9:
        return points
    axis = axis / dist
    
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    R = np.eye(3) + np.sin(angle_rad) * K + (1 - np.cos(angle_rad)) * (K @ K)
    return ((points - p1) @ R.T) + p1

def find_field(mesh, possible_names, is_cell=None):
    """ Finds a field in point or cell data and returns (data, associated_points) """
    if is_cell is True:
        data = mesh.cell_data
        pts = mesh.cell_centers().points
        for name in possible_names:
            if name in data.keys():
                return data[name], pts
    elif is_cell is False:
        data = mesh.point_data
        pts = mesh.points
        for name in possible_names:
            if name in data.keys():
                return data[name], pts
    else:
        # Search point data first
        for name in possible_names:
            if name in mesh.point_data.keys():
                return mesh.point_data[name], mesh.points
        # Then cell data
        for name in possible_names:
            if name in mesh.cell_data.keys():
                return mesh.cell_data[name], mesh.cell_centers().points
    return None, None

def extract_torso_landmarks(torso):
    tags, centers = find_field(torso, ['elemTag', 'tags', 'Material', 'material'], is_cell=True)
    if tags is None:
        raise ValueError(f"Torso missing tags. Point Data: {list(torso.point_data.keys())}, Cell Data: {list(torso.cell_data.keys())}")
    # 25=LV, 24=RV, 16-19=Valves
    lv_idx = np.where(tags == 25)[0]
    rv_idx = np.where(tags == 24)[0]
    
    if len(lv_idx) == 0 or len(rv_idx) == 0:
        print(f"  Warning: Torso tags 24/25 not found. Found: {np.unique(tags)}")

    lv_c = np.mean(centers[lv_idx], axis=0)
    rv_c = np.mean(centers[rv_idx], axis=0)
    
    valve_tags = [16, 17, 18, 19]
    valve_pts = [np.mean(centers[np.where(tags == vt)[0]], axis=0) 
                 for vt in valve_tags if len(np.where(tags == vt)[0]) > 0]
    base_c = np.mean(valve_pts, axis=0)
    
    vent_idx = np.where(np.isin(tags, [24, 25]))[0]
    vent_pts = centers[vent_idx]
    apex = vent_pts[np.argmin(vent_pts[:, 2])]
    
    landmarks = np.array([lv_c, rv_c, base_c, apex])
    print(f"  Torso Landmarks (LV, RV, Base, Apex):\n{landmarks}")
    return landmarks

def extract_heart_landmarks(heart):
    # Try multiple common naming patterns for robustness
    tags, tags_pts = find_field(heart, ['tags', 'Material', 'material', 'tag'])
    longi, longi_pts = find_field(heart, ['uvc_longitudinal', 'longitudinal', 'longi', 'z_coord'])
    
    if tags is None or longi is None:
        print("\n--- DEBUG: HEART MESH CONTENT ---")
        print(f"Point Data keys: {list(heart.point_data.keys())}")
        print(f"Cell Data keys: {list(heart.cell_data.keys())}")
        print("----------------------------------\n")
        raise ValueError("Heart mesh is missing 'tags' or 'uvc_longitudinal'. See DEBUG list above.")

    # User defined: -1=LV, 1=RV, 0=Apex, 1=Base
    lv_idx = np.where(tags == -1)[0]
    curr_tags_pts = tags_pts
    
    if len(lv_idx) == 0: 
        # Fallback to check if tags are 1/2 or if they are in infra
        infra, infra_pts = find_field(heart, ['uvc_intraventricular', 'intra', 'transmural'])
        if infra is not None:
            lv_idx = np.where(infra < -0.8)[0]
            rv_idx = np.where(infra > 0.8)[0]
            curr_tags_pts = infra_pts
        else:
            lv_idx = np.where(tags == 1)[0]
            rv_idx = np.where(tags == 2)[0]
            curr_tags_pts = tags_pts
    else:
        rv_idx = np.where(tags == 1)[0]

    lv_c = np.mean(curr_tags_pts[lv_idx], axis=0)
    rv_c = np.mean(curr_tags_pts[rv_idx], axis=0)
    
    base_idx = np.where(longi > 0.95)[0]
    apex_idx = np.where(longi < 0.05)[0]
    if len(base_idx) == 0: base_idx = [np.argmax(longi)]
    if len(apex_idx) == 0: apex_idx = [np.argmin(longi)]
    
    base_c = np.mean(longi_pts[base_idx], axis=0)
    apex_c = np.mean(longi_pts[apex_idx], axis=0)
    
    landmarks = np.array([lv_c, rv_c, base_c, apex_c])
    print(f"  Heart Landmarks (LV, RV, Base, Apex):\n{landmarks}")
    return landmarks

def main():
    torso_path = "KCL_torso1_scaled.vtu"
    heart_path = "ASCIIlegacy02_biventricular_conductivity.vtk"
    
    print(f"Loading {torso_path}...")
    torso = pv.read(torso_path)
    print(f"Loading {heart_path}...")
    heart = pv.read(heart_path)
    
    print("Extracting landmarks...")
    t_land = extract_torso_landmarks(torso)
    h_land = extract_heart_landmarks(heart)
    
    # APPLY LV/RV SWAP: Based on the "Switched" feedback
    # We map Heart[RV] to Torso[LV] and vice-versa to correct the inversion
    h_land_swapped = h_land.copy()
    h_land_swapped[0] = h_land[1] # Heart RV -> Torso LV position
    h_land_swapped[1] = h_land[0] # Heart LV -> Torso RV position
    
    print("Computing alignment with L/R swap...")
    R, t = rigid_transform_3D(h_land_swapped, t_land)
    
    print("Applying transformation...")
    pts = np.asarray(heart.points, dtype=np.float64)
    pts[~np.isfinite(pts)] = 0
    
    # --- CALCULATE TRANSFORMATIONS ---
    # manual parameters (from user feedback)
    manual_roll = 40      # degrees (roll adjustment)
    manual_tilt = -30     # degrees (tilt adjustment)
    
    if not np.all(np.isfinite(R)) or not np.all(np.isfinite(t)):
        raise ValueError("Rotation matrix or translation vector contains non-finite values.")
    
    print(f"Base Kabsch R:\n{R}")
    print(f"Base Kabsch t:\n{t}")
    
    # ---------------------------------------------------------
    # STRATEGY: MOVE TORSO TO HEART (INVERSE TRANSFORMATION)
    # This keeps the original heart mesh and seeds as the reference.
    # ---------------------------------------------------------
    
    print("\nTransforming torso to match original heart space...")
    pts_torso = np.asarray(torso.points, dtype=np.float64)
    
    # Pivot point in torso space for rotations
    base_torso = t_land[2]
    apex_torso = t_land[3]
    y_axis = np.array([0, 1, 0])
    
    # UNWIND THE TRANSFORMATIONS IN REVERSE ORDER:
    # 1. Undo Manual Tilt
    print(f"  Step 1: Undoing manual tilt ({-manual_tilt}deg)...")
    pts_torso = rotate_around_axis(pts_torso, base_torso, base_torso + y_axis, -manual_tilt)
    
    # 2. Undo 180deg Inversion + Manual Roll
    print(f"  Step 2: Undoing inversion + manual roll ({-(180 + manual_roll)}deg)...")
    pts_torso = rotate_around_axis(pts_torso, base_torso, apex_torso, -(180 + manual_roll))
    
    # 3. Undo Kabsch Alignment (P_heart = (P_torso - t) @ R)
    print("  Step 3: Undoing Kabsch alignment...")
    pts_torso = (pts_torso - t) @ R
    
    torso.points = pts_torso
    
    output_torso = "torso_aligned_to_heart.vtu"
    torso.save(output_torso)
    print(f"Success! Torso moved to heart space and saved to {output_torso}")
    
    # Optional: Transform Purkinje Graph if present (scaling only)
    purkinje_path = "biv-line-glued_m_graph.vtk"
    if Path(purkinje_path).exists():
        print(f"\nProcessing {purkinje_path}...")
        purkinje = pv.read(purkinje_path)
        # Scale by 1000 to match heart units (meter -> mm)
        print("  Scaling Purkinje graph by 1000...")
        purkinje.points = purkinje.points * 1000
        output_purkinje = "purkinje_scaled.vtk"
        purkinje.save(output_purkinje)
        print(f"  Scaled Purkinje graph saved to {output_purkinje}")

if __name__ == "__main__":
    main()
