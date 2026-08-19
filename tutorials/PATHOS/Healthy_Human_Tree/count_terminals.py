import argparse
import os
import sys
import numpy as np
import vtk
from analyze_aha_activation import load_case_endo, load_purkinje

def main():
    case_dir = "."
    time_dir = "0.05"
    radius = 5.0

    endo_pts, aha, lv = load_case_endo(case_dir)
    term_pts, act_times = load_purkinje(case_dir, time_dir)

    print(f"Total activated terminals loaded: {len(term_pts)}")

    lv_endo = endo_pts[lv]
    rv_endo = endo_pts[~lv]

    # LV Locator
    tp_lv = vtk.vtkPoints()
    for p in lv_endo:
        tp_lv.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
    poly_lv = vtk.vtkPolyData()
    poly_lv.SetPoints(tp_lv)
    loc_lv = vtk.vtkPointLocator()
    loc_lv.SetDataSet(poly_lv)
    loc_lv.BuildLocator()

    # RV Locator
    tp_rv = vtk.vtkPoints()
    for p in rv_endo:
        tp_rv.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
    poly_rv = vtk.vtkPolyData()
    poly_rv.SetPoints(tp_rv)
    loc_rv = vtk.vtkPointLocator()
    loc_rv.SetDataSet(poly_rv)
    loc_rv.BuildLocator()

    lv_count = 0
    rv_count = 0
    unassigned = 0

    for p in term_pts:
        j_lv = loc_lv.FindClosestPoint(float(p[0]), float(p[1]), float(p[2]))
        d_lv = ((p - lv_endo[j_lv]) ** 2).sum() ** 0.5
        
        j_rv = loc_rv.FindClosestPoint(float(p[0]), float(p[1]), float(p[2]))
        d_rv = ((p - rv_endo[j_rv]) ** 2).sum() ** 0.5

        if d_lv <= radius and d_lv <= d_rv:
            lv_count += 1
        elif d_rv <= radius and d_rv < d_lv:
            rv_count += 1
        else:
            unassigned += 1

    print(f"LV Terminals: {lv_count}")
    print(f"RV Terminals: {rv_count}")
    if unassigned > 0:
        print(f"Unassigned (distance > {radius}mm): {unassigned}")

if __name__ == "__main__":
    main()
