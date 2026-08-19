#!/usr/bin/env python3
"""
analyze_aha_activation.py

Calculates the average activation time of Purkinje terminals per AHA-17 segment.
Reads the Purkinje points and pvjNodes from `constant/purkinjeGraph`, and reads
the `activationTime` field from the specified time directory.

Usage:
    python3 analyze_aha_activation.py CASE_DIR TIME_DIR [--radius 5]
"""
import argparse
import os
import sys
import numpy as np

try:
    import vtk
except ImportError:
    sys.exit("Error: missing vtk module.")

import re
COUNT_RE = re.compile(r"\d+")
VECTOR_RE = re.compile(r"\(([^()]*)\)")
FACE_RE = re.compile(r"^\s*\d+\(([^()]*)\)")

def iter_foam_list(path, require_internal_field=False):
    seen_internal = not require_internal_field
    expected = None
    in_values = False
    index = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if not seen_internal:
                seen_internal = s.startswith("internalField") and "nonuniform" in s
                continue
            if not in_values:
                if expected is None and COUNT_RE.fullmatch(s):
                    expected = int(s)
                elif expected is not None and s == "(":
                    in_values = True
                continue
            if s == ")": return
            yield index, s
            index += 1
    sys.exit(f"Error: could not read OpenFOAM list from {path}.")

def read_label_list(path):
    return np.asarray([int(s) for _, s in iter_foam_list(path)], dtype=np.int64)

def read_points(path):
    points = []
    for _, s in iter_foam_list(path):
        m = VECTOR_RE.search(s)
        if m: points.append(np.fromstring(m.group(1), sep=" ", dtype=float))
    return np.asarray(points, dtype=float)

def read_selected_face_centres(path, points, face_ids):
    order = {int(fid): i for i, fid in enumerate(face_ids)}
    centres = np.empty((len(face_ids), 3), dtype=float)
    max_face = int(np.max(face_ids))
    for face_i, s in iter_foam_list(path):
        if face_i > max_face: break
        if face_i in order:
            m = FACE_RE.match(s)
            if m:
                point_ids = np.fromstring(m.group(1), sep=" ", dtype=np.int64)
                centres[order[face_i]] = points[point_ids].mean(axis=0)
    return centres

def read_selected_owners(path, face_ids):
    order = {int(fid): i for i, fid in enumerate(face_ids)}
    owners = np.empty(len(face_ids), dtype=np.int64)
    max_face = int(np.max(face_ids))
    for face_i, s in iter_foam_list(path):
        if face_i > max_face: break
        if face_i in order: owners[order[face_i]] = int(s)
    return owners

def read_selected_internal_scalars(path, indices, dtype=float):
    unique, inverse = np.unique(indices, return_inverse=True)
    selected = {int(cell_i): i for i, cell_i in enumerate(unique)}
    out = np.empty(len(unique), dtype=dtype)
    max_cell = int(unique[-1])
    for cell_i, s in iter_foam_list(path, require_internal_field=True):
        if cell_i > max_cell: break
        if cell_i in selected: out[selected[cell_i]] = dtype(s)
    return out[inverse]

def load_case_endo(case_dir):
    sets_dir = os.path.join(case_dir, "constant", "polyMesh", "sets")
    lv_path = os.path.join(sets_dir, "LVEndoFaces")
    rv_path = os.path.join(sets_dir, "RVEndoFaces")
    if not (os.path.exists(lv_path) and os.path.exists(rv_path)):
        sys.exit("Error: requires constant/polyMesh/sets/LVEndoFaces and RVEndoFaces.")

    lv_faces = read_label_list(lv_path)
    rv_faces = read_label_list(rv_path)
    face_ids = np.concatenate([lv_faces, rv_faces])

    mesh_dir = os.path.join(case_dir, "constant", "polyMesh")
    points = read_points(os.path.join(mesh_dir, "points"))
    pts = read_selected_face_centres(os.path.join(mesh_dir, "faces"), points, face_ids)
    owners = read_selected_owners(os.path.join(mesh_dir, "owner"), face_ids)

    aha_path = os.path.join(case_dir, "0", "AHA_Segment")
    if not os.path.exists(aha_path):
        sys.exit("Error: requires 0/AHA_Segment field.")
        
    aha = read_selected_internal_scalars(aha_path, owners, dtype=float).astype(int)

    lv = np.zeros(len(face_ids), dtype=bool)
    lv[:len(lv_faces)] = True
    return pts, aha, lv

def load_purkinje(case_dir, time_dir):
    graph_path = os.path.join(case_dir, "constant", "purkinjeGraph")
    if not os.path.exists(graph_path):
        sys.exit("Error: constant/purkinjeGraph not found.")
        
    # parse pvjNodes
    pvj_nodes = []
    points = []
    in_pvj = False
    in_points = False
    expected_pts = None
    with open(graph_path) as f:
        for line in f:
            s = line.strip()
            if s == "pvjNodes":
                in_pvj = True
            elif in_pvj and s == "(": continue
            elif in_pvj and s == ")": in_pvj = False
            elif in_pvj and COUNT_RE.fullmatch(s):
                pvj_nodes.append(int(s))
                
            elif s == "points":
                in_points = True
            elif in_points and expected_pts is None and COUNT_RE.fullmatch(s):
                expected_pts = int(s)
            elif in_points and s == "(": continue
            elif in_points and s == ")": in_points = False
            elif in_points:
                m = VECTOR_RE.search(s)
                if m: points.append(np.fromstring(m.group(1), sep=" ", dtype=float))
                
    points = np.asarray(points, dtype=float)
    
    act_path = os.path.join(case_dir, time_dir, "activationTime.purkinjeNetwork")
    if not os.path.exists(act_path):
        sys.exit(f"Error: {act_path} not found. Ensure simulation has run.")
        
    act_times = np.asarray([float(s) for _, s in iter_foam_list(act_path, require_internal_field=False)])
    
    term_pts = points[pvj_nodes]
    term_act = act_times[pvj_nodes]
    
    # filter unactivated points
    valid = term_act > 0.0
    return term_pts[valid], term_act[valid]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", help="OpenFOAM case directory")
    ap.add_argument("time", help="Time directory (e.g., 0.05)")
    ap.add_argument("--radius", type=float, default=5.0, help="Coverage radius in mm")
    args = ap.parse_args()

    endo_pts, aha, lv = load_case_endo(args.case)
    term_pts, act_times = load_purkinje(args.case, args.time)

    R = args.radius
    print(f"Loaded {len(term_pts)} activated terminals.")

    # Process LV AHA Segments (1-17)
    lv_endo = endo_pts[lv]
    lv_aha = aha[lv]
    
    tp = vtk.vtkPoints()
    for p in lv_endo:
        tp.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
    poly = vtk.vtkPolyData(); poly.SetPoints(tp)
    loc = vtk.vtkPointLocator(); loc.SetDataSet(poly); loc.BuildLocator()

    seg_act_times = {s: [] for s in range(1, 18)}
    
    for p, t in zip(term_pts, act_times):
        j = loc.FindClosestPoint(float(p[0]), float(p[1]), float(p[2]))
        d = ((p - lv_endo[j]) ** 2).sum() ** 0.5
        if d <= R:
            s = int(lv_aha[j])
            if 1 <= s <= 17:
                seg_act_times[s].append(t)

    print("\n=== Average Activation Time per LV AHA-17 Segment (ms) ===")
    def avg_t(s): 
        return f"{np.mean(seg_act_times[s])*1000:>6.1f}" if seg_act_times[s] else "   ---"
        
    print("                      Base      Mid     Apex")
    print(f"    Anterior:       {avg_t(1)}   {avg_t(7)}   {avg_t(13)}")
    print(f"    Anteroseptal:   {avg_t(2)}   {avg_t(8)}   {avg_t(14)}")
    print(f"    Inferoseptal:   {avg_t(3)}   {avg_t(9)}      ---")
    print(f"    Inferior:       {avg_t(4)}   {avg_t(10)}   {avg_t(15)}")
    print(f"    Inferolateral:  {avg_t(5)}   {avg_t(11)}   {avg_t(16)}")
    print(f"    Anterolateral:  {avg_t(6)}   {avg_t(12)}      ---")
    print(f"    Apex Cap:          ---      ---   {avg_t(17)}")

if __name__ == "__main__":
    main()
