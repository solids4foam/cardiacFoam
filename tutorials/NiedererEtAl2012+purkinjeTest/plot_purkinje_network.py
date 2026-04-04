"""
plot_purkinje_network.py
------------------------
Writes the idealised Purkinje Network defined in
    constant/electroProperties  (purkinjeNetwork subdict)
as a VTP (VTK XML PolyData) file that can be opened directly in ParaView.

The network is a bifurcating tree with 8 nodes and 7 edges injecting
activation into 4 PVJ (Purkinje-Ventricular Junction) sites spread
across the NiedererEtAl2012 slab mesh:

    root(0) ─── bifurc(1) ──── left-main(2) ──── left-ant(4)    @ ( 2.45, 1.15, 0.25) mm
                           │                 └─── left-post(5)   @ (12.45, 0.35, 0.65) mm
                           └─── right-main(3) ─── right-ant(6)  @ ( 7.45, 1.15, 0.25) mm
                                              └─── right-post(7) @ (17.45, 0.35, 0.65) mm

PVJ positions match the 3D mesh cell centroids supplied in pvjLocations.
The slab mesh (blockMeshDict) is 20 x 3 x 7 mm, 200x30x70 cells (dx=0.1 mm).

Output
------
    purkinjeNetwork.vtp   — lines + per-node / per-edge scalar data
    slabBoundingBox.vtp   — 12 edges showing the mesh outline (reference frame)

Usage
-----
    python plot_purkinje_network.py
    # then open purkinjeNetwork.vtp (and slabBoundingBox.vtp) in ParaView.
    # Colour nodes by 'nodeType': 0=root, 1=bifurcation, 2=terminal(PVJ)
    # Colour edges by 'edgeLength_mm'
"""

# ---------------------------------------------------------------------------
# Network definition  (mirrors constant/electroProperties exactly)
# ---------------------------------------------------------------------------

# Node coordinates [m] — placed to match their anatomical role in the slab.
# PVJ nodes (4-7) sit at the centroid of the corresponding 3D mesh cell.
#   cellID = i + j*200 + k*6000,  x=(i+0.5)*0.1e-3, y=(j+0.5)*0.1e-3, z=(k+0.5)*0.1e-3
#   node4: i=24,j=11,k=2  -> ( 2.45, 1.15, 0.25) mm  cellID=14224
#   node5: i=124,j=3,k=6  -> (12.45, 0.35, 0.65) mm  cellID=36724
#   node6: i=74,j=11,k=2  -> ( 7.45, 1.15, 0.25) mm  cellID=14274
#   node7: i=174,j=3,k=6  -> (17.45, 0.35, 0.65) mm  cellID=36774

NODE_COORDS = {          # (x, y, z) in metres
    0: ( 9.95e-3, 1.50e-3, 3.40e-3),   # root
    1: ( 9.95e-3, 1.50e-3, 2.60e-3),   # first bifurcation
    2: ( 7.45e-3, 1.50e-3, 1.60e-3),   # left  main branch point
    3: (12.45e-3, 1.50e-3, 1.60e-3),   # right main branch point
    4: ( 2.45e-3, 1.15e-3, 0.25e-3),   # PVJ
    5: (12.45e-3, 0.35e-3, 0.65e-3),   # PVJ
    6: ( 7.45e-3, 1.15e-3, 0.25e-3),   # PVJ
    7: (17.45e-3, 0.35e-3, 0.65e-3),   # PVJ
}

# Node type labels for colouring in ParaView
NODE_TYPE = {
    0: 0,   # root
    1: 1,   # bifurcation
    2: 1,   # bifurcation
    3: 1,   # bifurcation
    4: 2,   # PVJ terminal
    5: 2,   # PVJ terminal
    6: 2,   # PVJ terminal
    7: 2,   # PVJ terminal
}
NODE_TYPE_LEGEND = "0=root  1=bifurcation  2=PVJ-terminal"

NODE_LABELS = {
    0: "root",
    1: "bifurcation",
    2: "left-main",
    3: "right-main",
    4: "PVJ-left-ant",
    5: "PVJ-left-post",
    6: "PVJ-right-ant",
    7: "PVJ-right-post",
}

# Edges: (nodeA, nodeB, length_m, conductance_S_per_m)
EDGES = [
    (0, 1, 0.003, 3.0),   # root -> first bifurcation
    (1, 2, 0.005, 3.0),   # -> left  main
    (1, 3, 0.005, 3.0),   # -> right main
    (2, 4, 0.004, 3.0),   # left  anterior  terminal
    (2, 5, 0.004, 3.0),   # left  posterior terminal
    (3, 6, 0.004, 3.0),   # right anterior  terminal
    (3, 7, 0.004, 3.0),   # right posterior terminal
]

# Slab bounding box [m] (from blockMeshDict: scale 0.001, vertices (0,0,0)-(20,3,7) mm)
SLAB_MIN = (0.0,    0.0,    0.0)
SLAB_MAX = (20e-3,  3e-3,   7e-3)


# ---------------------------------------------------------------------------
# VTP writer helpers  (no external dependencies — pure stdlib)
# ---------------------------------------------------------------------------

def _flat(coords):
    """Flatten list of (x,y,z) tuples to a space-separated string."""
    return "\n          ".join(
        f"{x:.6e} {y:.6e} {z:.6e}" for (x, y, z) in coords
    )


def _ints(values):
    return " ".join(str(v) for v in values)


def _floats(values, fmt=".4f"):
    return " ".join(format(v, fmt) for v in values)


def write_network_vtp(path):
    """Write the Purkinje network as a VTP PolyData file."""
    n_nodes = len(NODE_COORDS)
    n_edges = len(EDGES)

    points  = [NODE_COORDS[i] for i in range(n_nodes)]
    conn    = []   # flattened connectivity: [a0,b0, a1,b1, ...]
    offsets = []   # cumulative end-of-line indices
    edge_lengths_mm   = []
    edge_conductances = []

    for i, (a, b, L, sigma) in enumerate(EDGES):
        conn.extend([a, b])
        offsets.append(2 * (i + 1))
        edge_lengths_mm.append(L * 1e3)
        edge_conductances.append(sigma)

    node_types = [NODE_TYPE[i] for i in range(n_nodes)]

    xml = f"""\
<?xml version="1.0"?>
<!-- Purkinje Network — NiedererEtAl2012+purkinjeTest
     Slab mesh: 20x3x7 mm, 200x30x70 cells (dx=0.1 mm)
     Node types: {NODE_TYPE_LEGEND}
     Units: coordinates [m], edgeLength [mm], conductance [S/m]
-->
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{n_nodes}" NumberOfVerts="0"
           NumberOfLines="{n_edges}" NumberOfStrips="0" NumberOfPolys="0">

      <!-- ===== Node (point) coordinates ===== -->
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          {_flat(points)}
        </DataArray>
      </Points>

      <!-- ===== Edge (line) connectivity ===== -->
      <Lines>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {_ints(conn)}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {_ints(offsets)}
        </DataArray>
      </Lines>

      <!-- ===== Per-node data ===== -->
      <PointData Scalars="nodeType">
        <!-- nodeType: 0=root  1=bifurcation  2=PVJ-terminal -->
        <DataArray type="Int32" Name="nodeType" format="ascii">
          {_ints(node_types)}
        </DataArray>
      </PointData>

      <!-- ===== Per-edge data ===== -->
      <CellData Scalars="edgeLength_mm">
        <DataArray type="Float64" Name="edgeLength_mm" format="ascii">
          {_floats(edge_lengths_mm)}
        </DataArray>
        <DataArray type="Float64" Name="conductance_S_per_m" format="ascii">
          {_floats(edge_conductances)}
        </DataArray>
      </CellData>

    </Piece>
  </PolyData>
</VTKFile>
"""
    with open(path, "w") as f:
        f.write(xml)
    print(f"Written: {path}")


def write_slab_vtp(path):
    """Write the mesh bounding box as 12 line edges (reference frame)."""
    xmin, ymin, zmin = SLAB_MIN
    xmax, ymax, zmax = SLAB_MAX

    # 8 corners of the box
    corners = [
        (xmin, ymin, zmin), (xmax, ymin, zmin),
        (xmax, ymax, zmin), (xmin, ymax, zmin),
        (xmin, ymin, zmax), (xmax, ymin, zmax),
        (xmax, ymax, zmax), (xmin, ymax, zmax),
    ]

    # 12 edges of the box
    box_edges = [
        (0,1),(1,2),(2,3),(3,0),   # bottom face
        (4,5),(5,6),(6,7),(7,4),   # top face
        (0,4),(1,5),(2,6),(3,7),   # vertical edges
    ]

    conn    = []
    offsets = []
    for i, (a, b) in enumerate(box_edges):
        conn.extend([a, b])
        offsets.append(2 * (i + 1))

    xml = f"""\
<?xml version="1.0"?>
<!-- Slab bounding box: {int(xmax*1e3)}x{int(ymax*1e3)}x{int(zmax*1e3)} mm
     Open in ParaView alongside purkinjeNetwork.vtp for reference frame. -->
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="8" NumberOfVerts="0"
           NumberOfLines="12" NumberOfStrips="0" NumberOfPolys="0">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          {_flat(corners)}
        </DataArray>
      </Points>
      <Lines>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {_ints(conn)}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {_ints(offsets)}
        </DataArray>
      </Lines>
    </Piece>
  </PolyData>
</VTKFile>
"""
    with open(path, "w") as f:
        f.write(xml)
    print(f"Written: {path}")


# ---------------------------------------------------------------------------
# Print a human-readable summary to stdout
# ---------------------------------------------------------------------------

def print_summary():
    print("\n=== Purkinje Network Summary ===")
    print(f"  Nodes : {len(NODE_COORDS)}")
    print(f"  Edges : {len(EDGES)}")
    print(f"  PVJs  : {sum(1 for t in NODE_TYPE.values() if t == 2)}\n")
    print(f"  {'Node':<5}  {'Label':<20}  {'x (mm)':>8}  {'y (mm)':>8}  {'z (mm)':>8}  {'Type'}")
    print(f"  {'-'*5}  {'-'*20}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*15}")
    for i, (x, y, z) in NODE_COORDS.items():
        tname = ["root", "bifurcation", "PVJ-terminal"][NODE_TYPE[i]]
        print(f"  {i:<5}  {NODE_LABELS[i]:<20}  {x*1e3:8.2f}  {y*1e3:8.2f}  {z*1e3:8.2f}  {tname}")
    print()
    print(f"  {'Edge':<5}  {'A->B':<6}  {'L (mm)':>8}  {'sigma (S/m)':>12}")
    print(f"  {'-'*5}  {'-'*6}  {'-'*8}  {'-'*12}")
    for i, (a, b, L, sigma) in enumerate(EDGES):
        print(f"  {i:<5}  {a}->{b:<3}  {L*1e3:8.2f}  {sigma:12.4f}")
    print()
    print(f"  Slab bounding box: "
          f"X=[{SLAB_MIN[0]*1e3:.0f}, {SLAB_MAX[0]*1e3:.0f}] mm  "
          f"Y=[{SLAB_MIN[1]*1e3:.0f}, {SLAB_MAX[1]*1e3:.0f}] mm  "
          f"Z=[{SLAB_MIN[2]*1e3:.0f}, {SLAB_MAX[2]*1e3:.0f}] mm")
    print("================================\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os
    out_dir = os.path.dirname(os.path.abspath(__file__))
    constant_dir = os.path.join(out_dir, "constant")
    print_summary()
    write_network_vtp(os.path.join(constant_dir, "purkinjeNetwork.vtp"))
    write_slab_vtp(os.path.join(out_dir, "slabBoundingBox.vtp"))
    print("\nOpen both files in ParaView.")
    print("Suggested display:")
    print("  constant/purkinjeNetwork.vtp  -> Tube filter, colour by 'nodeType'")
    print("  slabBoundingBox.vtp  -> Wireframe, opacity 0.3")
