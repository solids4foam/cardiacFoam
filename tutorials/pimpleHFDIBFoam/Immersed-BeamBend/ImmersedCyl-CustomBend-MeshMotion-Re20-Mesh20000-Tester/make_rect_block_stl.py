import os
import math

# --------------------------------------------------
# Geometry settings
# --------------------------------------------------
Lx = 24.0
Ly = 4.1

# Beam dimensions
block_x = 0.2                 # thickness in x
block_y = 0.495 * Ly          # height in y
z0 = -0.05
z1 =  0.05                    # thickness in z = 0.1

# Beam location (centered at x = 12)
x_center = 12.0
x0 = x_center - 0.5 * block_x
x1 = x_center + 0.5 * block_x
y0 = 0.0
y1 = block_y

# Number of segments along beam height
Ny = 40

outdir = os.path.join("constant", "triSurface")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "rectBlock_centered.stl")

# --------------------------------------------------
# Helpers
# --------------------------------------------------
def normal(p1, p2, p3):
    ux, uy, uz = p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]
    vx, vy, vz = p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    mag = math.sqrt(nx*nx + ny*ny + nz*nz)
    if mag == 0.0:
        return (0.0, 0.0, 0.0)
    return (nx/mag, ny/mag, nz/mag)

def write_facet(f, p1, p2, p3):
    n = normal(p1, p2, p3)
    f.write(f"  facet normal {n[0]:.8e} {n[1]:.8e} {n[2]:.8e}\n")
    f.write("    outer loop\n")
    f.write(f"      vertex {p1[0]:.8e} {p1[1]:.8e} {p1[2]:.8e}\n")
    f.write(f"      vertex {p2[0]:.8e} {p2[1]:.8e} {p2[2]:.8e}\n")
    f.write(f"      vertex {p3[0]:.8e} {p3[1]:.8e} {p3[2]:.8e}\n")
    f.write("    endloop\n")
    f.write("  endfacet\n")

# --------------------------------------------------
# Build segmented beam
# --------------------------------------------------
ys = [y0 + (y1 - y0) * i / Ny for i in range(Ny + 1)]

# Points on the four vertical edges
left_back  = [(x0, y, z0) for y in ys]
right_back = [(x1, y, z0) for y in ys]
left_front = [(x0, y, z1) for y in ys]
right_front= [(x1, y, z1) for y in ys]

with open(outfile, "w") as f:
    f.write("solid rectBlock_centered\n")

    # ------------------------------
    # Back face (z = z0), outward = -z
    # ------------------------------
    for i in range(Ny):
        p00 = left_back[i]
        p10 = right_back[i]
        p11 = right_back[i+1]
        p01 = left_back[i+1]

        write_facet(f, p00, p11, p10)
        write_facet(f, p00, p01, p11)

    # ------------------------------
    # Front face (z = z1), outward = +z
    # ------------------------------
    for i in range(Ny):
        p00 = left_front[i]
        p10 = right_front[i]
        p11 = right_front[i+1]
        p01 = left_front[i+1]

        write_facet(f, p00, p10, p11)
        write_facet(f, p00, p11, p01)

    # ------------------------------
    # Left face (x = x0), outward = -x
    # ------------------------------
    for i in range(Ny):
        p00 = left_back[i]
        p10 = left_front[i]
        p11 = left_front[i+1]
        p01 = left_back[i+1]

        write_facet(f, p00, p10, p11)
        write_facet(f, p00, p11, p01)

    # ------------------------------
    # Right face (x = x1), outward = +x
    # ------------------------------
    for i in range(Ny):
        p00 = right_back[i]
        p10 = right_front[i]
        p11 = right_front[i+1]
        p01 = right_back[i+1]

        write_facet(f, p00, p11, p10)
        write_facet(f, p00, p01, p11)

    # ------------------------------
    # Bottom face (y = y0), outward = -y
    # ------------------------------
    p0 = (x0, y0, z0)
    p1 = (x1, y0, z0)
    p2 = (x1, y0, z1)
    p3 = (x0, y0, z1)
    write_facet(f, p0, p2, p1)
    write_facet(f, p0, p3, p2)

    # ------------------------------
    # Top face (y = y1), outward = +y
    # ------------------------------
    p0 = (x0, y1, z0)
    p1 = (x1, y1, z0)
    p2 = (x1, y1, z1)
    p3 = (x0, y1, z1)
    write_facet(f, p0, p1, p2)
    write_facet(f, p0, p2, p3)

    f.write("endsolid rectBlock_centered\n")

print("STL written to:", outfile)
print(f"x = [{x0}, {x1}]  thickness = {x1 - x0}")
print(f"y = [{y0}, {y1}]  height = {y1 - y0}")
print(f"z = [{z0}, {z1}]  thickness = {z1 - z0}")
print(f"Ny segments = {Ny}")
print(f"Expected volume = {(x1-x0)*(y1-y0)*(z1-z0):.8f}")
print(f"Expected triangles = {8*Ny + 4}")
