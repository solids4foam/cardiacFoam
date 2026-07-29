import sys
import subprocess
import re
import argparse
from pathlib import Path

def characterize_mesh(case_dir):
    case_path = Path(case_dir)
    if not case_path.exists():
        print(f"Error: {case_dir} does not exist.")
        return None
    
    print(f"Running checkMesh on {case_dir}...")
    
    # Run checkMesh. We assume the environment is sourced.
    try:
        result = subprocess.run(
            ["checkMesh", "-time", "0"], 
            cwd=case_path, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
    except FileNotFoundError:
        print("Error: checkMesh command not found. Ensure OpenFOAM environment is sourced.")
        return None
        
    output = result.stdout
    
    # Extract metrics using regex
    metrics = {
        "Name": case_path.name,
        "Points": None,
        "Cells": None,
        "Faces": None,
        "Internal Faces": None,
        "Bounding Box": None,
        "Max Non-Orthogonality": None,
        "Average Non-Orthogonality": None,
        "Severely Non-Orthogonal Faces": 0,
        "Max Skewness": None,
        "Highly Skew Faces": 0,
        "Tetrahedra": 0,
        "Polyhedra": 0
    }
    
    # Parsers
    m_points = re.search(r"points:\s+(\d+)", output)
    if m_points: metrics["Points"] = int(m_points.group(1))
        
    m_cells = re.search(r"cells:\s+(\d+)", output)
    if m_cells: metrics["Cells"] = int(m_cells.group(1))
        
    m_faces = re.search(r"faces:\s+(\d+)", output)
    if m_faces: metrics["Faces"] = int(m_faces.group(1))
        
    m_int_faces = re.search(r"internal faces:\s+(\d+)", output)
    if m_int_faces: metrics["Internal Faces"] = int(m_int_faces.group(1))
        
    m_bbox = re.search(r"Overall domain bounding box\s+(.*)", output)
    if m_bbox: metrics["Bounding Box"] = m_bbox.group(1).strip()
        
    m_nonortho = re.search(r"Mesh non-orthogonality Max:\s+([0-9.]+)\s+average:\s+([0-9.]+)", output)
    if m_nonortho:
        metrics["Max Non-Orthogonality"] = float(m_nonortho.group(1))
        metrics["Average Non-Orthogonality"] = float(m_nonortho.group(2))
        
    m_sev_nonortho = re.search(r"severely non-orthogonal \(> \d+ degrees\) faces:\s+(\d+)", output)
    if m_sev_nonortho: metrics["Severely Non-Orthogonal Faces"] = int(m_sev_nonortho.group(1))
        
    m_skew = re.search(r"Max skewness = ([0-9.]+).*?(\d+) highly skew faces", output)
    if m_skew:
        metrics["Max Skewness"] = float(m_skew.group(1))
        metrics["Highly Skew Faces"] = int(m_skew.group(2))
        
    m_tet = re.search(r"tetrahedra:\s+(\d+)", output)
    if m_tet: metrics["Tetrahedra"] = int(m_tet.group(1))
        
    m_poly = re.search(r"polyhedra:\s+(\d+)", output)
    if m_poly: metrics["Polyhedra"] = int(m_poly.group(1))
        
    return metrics

def main():
    parser = argparse.ArgumentParser(description="Characterize OpenFOAM meshes")
    parser.add_argument("--cases", nargs="+", default=["monodomain_coarse", "monodomain_fine"], help="Case directories to check")
    args = parser.parse_args()
    
    results = []
    for case in args.cases:
        res = characterize_mesh(case)
        if res:
            results.append(res)
            
    if not results:
        print("No valid results found.")
        return
        
    # Print Markdown Table
    print("\n# Mesh Characterization Report\n")
    headers = [
        "Case", "Cells (Polyhedra)", "Points", "Faces", 
        "Max Non-Ortho", "Avg Non-Ortho", "Severe Non-Ortho",
        "Max Skewness", "Highly Skew"
    ]
    
    # Table formatting
    print("| " + " | ".join(headers) + " |")
    print("|" + "|".join(["---"] * len(headers)) + "|")
    
    for r in results:
        cells_str = f"{r['Cells']} ({r['Polyhedra']})"
        row = [
            r["Name"],
            cells_str,
            str(r["Points"]),
            str(r["Faces"]),
            f"{r['Max Non-Orthogonality']:.2f}" if r['Max Non-Orthogonality'] else "N/A",
            f"{r['Average Non-Orthogonality']:.2f}" if r['Average Non-Orthogonality'] else "N/A",
            str(r["Severely Non-Orthogonal Faces"]),
            f"{r['Max Skewness']:.2f}" if r['Max Skewness'] else "N/A",
            str(r["Highly Skew Faces"])
        ]
        print("| " + " | ".join(row) + " |")

if __name__ == "__main__":
    main()
