import math

def write_ascii_stl(filename, triangles, solid_name="cylinder"):
    with open(filename, "w") as f:
        f.write(f"solid {solid_name}\n")
        for (nx, ny, nz), (v1, v2, v3) in triangles:
            f.write(f"  facet normal {nx:.8e} {ny:.8e} {nz:.8e}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {v1[0]:.8e} {v1[1]:.8e} {v1[2]:.8e}\n")
            f.write(f"      vertex {v2[0]:.8e} {v2[1]:.8e} {v2[2]:.8e}\n")
            f.write(f"      vertex {v3[0]:.8e} {v3[1]:.8e} {v3[2]:.8e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")

def normal_of_triangle(v1, v2, v3):
    ax, ay, az = (v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2])
    bx, by, bz = (v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2])
    nx = ay*bz - az*by
    ny = az*bx - ax*bz
    nz = ax*by - ay*bx
    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
    if norm < 1e-16:
        return (0.0, 0.0, 0.0)
    return (nx/norm, ny/norm, nz/norm)

def make_cylinder(radius=0.5, zmin=-0.261725, zmax=0.261725, ntheta=128):
    triangles = []

    # Ring points
    thetas = [2*math.pi*i/ntheta for i in range(ntheta)]
    ring_bot = [(radius*math.cos(t), radius*math.sin(t), zmin) for t in thetas]
    ring_top = [(radius*math.cos(t), radius*math.sin(t), zmax) for t in thetas]

    # Side surface (two triangles per segment)
    for i in range(ntheta):
        j = (i + 1) % ntheta
        v1 = ring_bot[i]
        v2 = ring_bot[j]
        v3 = ring_top[j]
        v4 = ring_top[i]
        # triangle (v1,v2,v3)
        n = normal_of_triangle(v1, v2, v3)
        triangles.append((n, (v1, v2, v3)))
        # triangle (v1,v3,v4)
        n = normal_of_triangle(v1, v3, v4)
        triangles.append((n, (v1, v3, v4)))

    # Caps (fan triangulation)
    c_bot = (0.0, 0.0, zmin)
    c_top = (0.0, 0.0, zmax)

    # Bottom cap normal should point -z (outside)
    for i in range(ntheta):
        j = (i + 1) % ntheta
        v1 = c_bot
        v2 = ring_bot[j]
        v3 = ring_bot[i]
        n = normal_of_triangle(v1, v2, v3)
        triangles.append((n, (v1, v2, v3)))

    # Top cap normal should point +z (outside)
    for i in range(ntheta):
        j = (i + 1) % ntheta
        v1 = c_top
        v2 = ring_top[i]
        v3 = ring_top[j]
        n = normal_of_triangle(v1, v2, v3)
        triangles.append((n, (v1, v2, v3)))

    return triangles

if __name__ == "__main__":
    tris = make_cylinder(radius=0.5, zmin=-0.261725, zmax=0.261725, ntheta=128)
    out = "constant/triSurface/cylinder_R0p5_L0p52345_center0.stl"
    write_ascii_stl(out, tris, solid_name="cylinder")
    print(f"Wrote STL: {out}  (triangles: {len(tris)})")
