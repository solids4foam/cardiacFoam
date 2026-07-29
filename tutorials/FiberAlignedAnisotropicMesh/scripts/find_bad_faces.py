import sys
import numpy as np
import re
from pathlib import Path

def parse_openfoam_set(filepath):
    """Parses an OpenFOAM set file and returns a list of integer IDs."""
    ids = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        in_list = False
        for line in lines:
            line = line.strip()
            if not line or line.startswith('//') or line.startswith('/*'):
                continue
            if line.isdigit() and len(ids) == 0:
                continue # size of the list
            if line == '(':
                in_list = True
                continue
            if line == ')' or line == ');':
                in_list = False
                continue
            if in_list and line.isdigit():
                ids.append(int(line))
    except Exception as e:
        pass
    return ids

def parse_openfoam_faces(case_dir):
    """Parses the faces file and points file to compute face centers."""
    points_file = Path(case_dir) / "constant" / "polyMesh" / "points"
    faces_file = Path(case_dir) / "constant" / "polyMesh" / "faces"
    
    # Very basic parser for points
    points = []
    try:
        with open(points_file, 'r') as f:
            content = f.read()
            # Extract the points array
            m = re.search(r'\(\s*(.*?)\s*\)', content, re.DOTALL)
            if m:
                # Find all (x y z)
                pts_str = re.findall(r'\((.*?)\)', m.group(1))
                for pt in pts_str:
                    coords = list(map(float, pt.split()))
                    points.append(coords)
    except Exception:
        pass
        
    # Basic parser for faces
    faces = []
    try:
        with open(faces_file, 'r') as f:
            content = f.read()
            m = re.search(r'\(\s*(.*?)\s*\)', content, re.DOTALL)
            if m:
                face_lines = m.group(1).strip().split('\n')
                for line in face_lines:
                    line = line.strip()
                    if line.startswith('3(') or line.startswith('4(') or line.startswith('5(') or ')' in line:
                        # Extract the numbers inside the parenthesis
                        fm = re.search(r'\((.*?)\)', line)
                        if fm:
                            p_ids = list(map(int, fm.group(1).split()))
                            faces.append(p_ids)
    except Exception:
        pass
        
    return points, faces

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 find_bad_faces.py <case_dir>")
        sys.exit(1)
        
    case_dir = Path(sys.argv[1])
    
    # 1. Run checkMesh first
    print(f"Checking {case_dir.name}...")
    
    # 2. Parse bad faces sets
    non_ortho = parse_openfoam_set(case_dir / "constant" / "polyMesh" / "sets" / "nonOrthoFaces")
    skew = parse_openfoam_set(case_dir / "constant" / "polyMesh" / "sets" / "skewFaces")
    wrong_orient = parse_openfoam_set(case_dir / "constant" / "polyMesh" / "sets" / "wrongOrientedFaces")
    
    if not non_ortho and not skew and not wrong_orient:
        print("No bad face sets found. Did you run checkMesh on this case?")
        return
        
    print(f"Found {len(non_ortho)} non-orthogonal faces.")
    print(f"Found {len(skew)} highly skew faces.")
    print(f"Found {len(wrong_orient)} wrongly oriented faces.")
    
    # 3. Try to parse geometry to get coordinates
    points, faces = parse_openfoam_faces(case_dir)
    if not points or not faces:
        print("Could not parse points/faces to find geometric locations.")
        print("Use: 'foamToVTK -faceSet nonOrthoFaces' instead to view them in ParaView.")
        return
        
    def print_locations(face_ids, name):
        if not face_ids: return
        print(f"\n--- {name} Locations (X, Y, Z) ---")
        for fid in face_ids:
            if fid < len(faces):
                f_pts = faces[fid]
                coords = [points[pid] for pid in f_pts if pid < len(points)]
                if coords:
                    center = np.mean(coords, axis=0)
                    print(f"Face {fid}: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
                    
    print_locations(non_ortho[:10], "Non-Orthogonal Faces (First 10)")
    print_locations(skew[:10], "Skew Faces (First 10)")
    print_locations(wrong_orient[:10], "Wrongly Oriented Faces (First 10)")
    
    print("\nTip: Run 'foamToVTK -faceSet nonOrthoFaces' in the case directory to visualize these in ParaView.")

if __name__ == "__main__":
    main()
