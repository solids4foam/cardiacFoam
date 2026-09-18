#!/usr/bin/env python3
import math

def cross(ax, ay, az, bx, by, bz):
    return (
        ay*bz - az*by,
        az*bx - ax*bz,
        ax*by - ay*bx
    )

def norm(nx, ny, nz):
    mag = math.sqrt(nx*nx + ny*ny + nz*nz)
    if mag == 0.0:
        return (0.0, 0.0, 0.0)
    return (nx/mag, ny/mag, nz/mag)

def facet(f, v1, v2, v3):
    # Compute normal from (v2-v1) x (v3-v1)
    ax, ay, az = (v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2])
    bx, by, bz = (v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2])
    nx, ny, nz = cross(ax, ay, az, bx, by, bz)
    nx, ny, nz = norm(nx, ny, nz)

    f.write(f"  facet normal {nx:.8e} {ny:.8e} {nz:.8e}\n")
    f.write("    outer loop\n")
    f.write(f"      vertex {v1[0]:.8e} {v1[1]:.8e} {v1[2]:.8e}\n")
    f.write(f"      vertex {v2[0]:.8e} {v2[1]:.8e} {v2[2]:.8e}\n")
    f.write(f"      vertex {v3[0]:.8e} {v3[1]:.8e} {v3[2]:.8e}\n")
    f.write("    endloop\n")
    f.write("  endfacet\n")

def write_cylinder_stl(
    filename="cylinder_fits_old_hole.stl",
    cx=1.1, cy=0.2,
    radius=0.05,
    zmin=0.0, zmax=0.1,
    nTheta=180,
    cap_bottom=True,
    cap_top=True
):
    # Precompute circle points
    ptsB = []
    ptsT = []
    for i in range(nTheta):
        th = 2.0*math.pi*i/nTheta
        x = cx + radius*math.cos(th)
        y = cy + radius*math.sin(th)
        ptsB.append((x, y, zmin))
        ptsT.append((x, y, zmax))

    with open(filename, "w") as f:
        f.write("solid cylinder\n")

        # Side surface: 2 triangles per segment
        for i in range(nTheta):
            j = (i + 1) % nTheta
            v_ib = ptsB[i]
            v_jb = ptsB[j]
            v_it = ptsT[i]
            v_jt = ptsT[j]

            # Outward normals with this winding (right-hand rule)
            facet(f, v_ib, v_jb, v_jt)
            facet(f, v_ib, v_jt, v_it)

        # Bottom cap (normal outward is -z): reverse winding
        if cap_bottom:
            cB = (cx, cy, zmin)
            for i in range(nTheta):
                j = (i + 1) % nTheta
                facet(f, cB, ptsB[j], ptsB[i])

        # Top cap (normal outward is +z)
        if cap_top:
            cT = (cx, cy, zmax)
            for i in range(nTheta):
                j = (i + 1) % nTheta
                facet(f, cT, ptsT[i], ptsT[j])

        f.write("endsolid cylinder\n")

    print(f"Wrote {filename}")
    print(f"centre=({cx}, {cy}, {(zmin+zmax)/2.0}), radius={radius}, z=[{zmin},{zmax}], nTheta={nTheta}")

if __name__ == "__main__":
    # EXACT match to your old hole location/size:
    write_cylinder_stl(
        filename="cylinder_fits_old_hole.stl",
        cx=1.1, cy=0.2,
        radius=0.05,
        zmin=0.0, zmax=0.1,
        nTheta=180,
        cap_bottom=True,
        cap_top=True
    )
