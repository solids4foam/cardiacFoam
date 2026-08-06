import sys

def scale_vtk(in_path, out_path, scale):
    with open(in_path, 'r') as f:
        lines = f.readlines()
        
    out = []
    in_points = False
    points_to_read = 0
    
    for line in lines:
        if line.startswith("POINTS"):
            in_points = True
            parts = line.strip().split()
            points_to_read = int(parts[1])
            out.append(line)
            continue
            
        if in_points and points_to_read > 0:
            parts = line.strip().split()
            if len(parts) >= 3:
                # Some lines might have multiple points? Usually 3 floats per point per line.
                # Just scale all floats on the line
                scaled_floats = [str(float(x) * scale) for x in parts]
                out.append(" ".join(scaled_floats) + "\n")
                # usually VTK ascii from openfoam is 1 point per line or 3 floats per line
                # wait, if there are multiple points per line, this scales them all
                points_to_read -= len(parts) // 3
            else:
                out.append(line)
            if points_to_read <= 0:
                in_points = False
            continue
            
        out.append(line)
        
    with open(out_path, 'w') as f:
        f.writelines(out)

if __name__ == "__main__":
    scale_vtk("purkinje.vtk", "purkinje_scaled.vtk", 0.001)
    print("Manual ASCII scale complete.")
