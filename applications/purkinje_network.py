#!/usr/bin/env python3
"""Generate Purkinje input files for OpenFOAM.

Primary output:
    purkinje_network.msh
    purkinje_graph.vtk

Use this with:
    gmshToFoam purkinje_network.msh
    1DgraphToFoam purkinje_graph.vtk

The .msh file is the optional tagged tube mesh. The line-only graph VTK is the
input consumed by 1DgraphToFoam for the graph-based conductionSystemDomain.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


RADIUS = 0.0002
MESH_SIZE = 0.0003

TAG_ROOT = 1
TAG_TERMINALS = 2
TAG_WALL = 3
TAG_VOLUME = 4


NODES = [
    [0.001, 0.0015, 0.0035],  # 0: His bundle entry (root)
    [0.005, 0.0015, 0.0035],  # 1: Trunk bifurcation
    [0.008, 0.0015, 0.0055],  # 2: Left fascicle
    [0.008, 0.0015, 0.0015],  # 3: Right fascicle
    [0.011, 0.0015, 0.0060],  # 4: Sub-branch upper-left
    [0.011, 0.0015, 0.0050],  # 5: Sub-branch upper-right
    [0.011, 0.0015, 0.0020],  # 6: Sub-branch lower-left
    [0.011, 0.0015, 0.0010],  # 7: Sub-branch lower-right
    [0.014, 0.0015, 0.0060],  # 8: Terminal
    [0.014, 0.0015, 0.0055],  # 9: Terminal
    [0.014, 0.0015, 0.0050],  # 10: Terminal
    [0.014, 0.0015, 0.0040],  # 11: Terminal
    [0.014, 0.0015, 0.0030],  # 12: Terminal
    [0.014, 0.0015, 0.0020],  # 13: Terminal
    [0.017, 0.0015, 0.0015],  # 14: Terminal
    [0.017, 0.0015, 0.0010],  # 15: Terminal
]


EDGES = [
    [0, 1],
    [1, 2],
    [2, 4],
    [2, 5],
    [1, 3],
    [3, 6],
    [3, 7],
    [4, 8],
    [4, 9],
    [5, 10],
    [5, 11],
    [6, 12],
    [6, 13],
    [7, 14],
    [7, 15],
]


TERMINAL_NODE_IDS = tuple(range(8, 16))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a tagged 3D Purkinje network mesh."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where purkinje_network.msh and .vtk are written.",
    )
    parser.add_argument(
        "--basename",
        default="purkinje_network",
        help="Output basename without extension.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=RADIUS,
        help="Purkinje tube radius in meters.",
    )
    parser.add_argument(
        "--mesh-size",
        type=float,
        default=MESH_SIZE,
        help="Target 3D mesh element size in meters.",
    )
    parser.add_argument(
        "--graph-step",
        type=float,
        default=MESH_SIZE,
        help="Maximum segment length in the line graph VTK.",
    )
    parser.add_argument(
        "--diffusivity",
        type=float,
        default=1.0,
        help="Uniform line diffusivity/conductance value written as CELL_DATA.",
    )
    parser.add_argument(
        "--graph-only",
        action="store_true",
        help="Only write the line graph VTK; skip gmsh and the 3D tube mesh.",
    )
    return parser.parse_args()


def add_tube_volume(nodes: list[list[float]], edges: list[list[int]], radius: float):
    volumes = []

    for point in nodes:
        tag = gmsh.model.occ.addSphere(point[0], point[1], point[2], radius)
        volumes.append((3, tag))

    for start_i, end_i in edges:
        start = nodes[start_i]
        end = nodes[end_i]
        direction = [
            end[0] - start[0],
            end[1] - start[1],
            end[2] - start[2],
        ]
        tag = gmsh.model.occ.addCylinder(
            start[0],
            start[1],
            start[2],
            direction[0],
            direction[1],
            direction[2],
            radius,
        )
        volumes.append((3, tag))

    fused, _ = gmsh.model.occ.fuse([volumes[0]], volumes[1:])
    gmsh.model.occ.synchronize()
    return fused


def surface_tags_near_point(point: list[float], radius: float) -> list[int]:
    eps = 1.1 * radius
    surfaces = gmsh.model.getEntitiesInBoundingBox(
        point[0] - eps,
        point[1] - eps,
        point[2] - eps,
        point[0] + eps,
        point[1] + eps,
        point[2] + eps,
        dim=2,
    )
    return [tag for _, tag in surfaces]


def add_physical_groups(nodes: list[list[float]], radius: float) -> None:
    root_tags = surface_tags_near_point(nodes[0], radius)

    terminal_tags: list[int] = []
    for node_i in TERMINAL_NODE_IDS:
        terminal_tags.extend(surface_tags_near_point(nodes[node_i], radius))
    terminal_tags = sorted(set(terminal_tags) - set(root_tags))

    all_surface_tags = [tag for _, tag in gmsh.model.getEntities(dim=2)]
    wall_tags = sorted(set(all_surface_tags) - set(root_tags) - set(terminal_tags))

    volume_tags = [tag for _, tag in gmsh.model.getEntities(dim=3)]

    required_groups = {
        "root": root_tags,
        "terminals": terminal_tags,
        "purkinje_wall": wall_tags,
        "purkinje_volume": volume_tags,
    }
    missing = [name for name, tags in required_groups.items() if not tags]
    if missing:
        raise RuntimeError(f"Empty physical group(s): {', '.join(missing)}")

    gmsh.model.addPhysicalGroup(2, root_tags, tag=TAG_ROOT, name="root")
    gmsh.model.addPhysicalGroup(
        2,
        terminal_tags,
        tag=TAG_TERMINALS,
        name="terminals",
    )
    gmsh.model.addPhysicalGroup(2, wall_tags, tag=TAG_WALL, name="purkinje_wall")
    gmsh.model.addPhysicalGroup(3, volume_tags, tag=TAG_VOLUME, name="purkinje_volume")

    print("Physical groups:")
    for name, tags in required_groups.items():
        print(f"  {name}: {len(tags)} entity/entities")


def write_graph_vtk(
    output_dir: Path,
    basename: str,
    nodes: list[list[float]],
    edges: list[list[int]],
    terminal_node_ids: tuple[int, ...],
    radius: float,
    graph_step: float,
    diffusivity: float,
) -> None:
    graph_path = output_dir / f"{basename}_graph.vtk"

    graph_points = [tuple(point) for point in nodes]
    node_roles = [0 for _ in graph_points]
    node_roles[0] = 1
    for node_i in terminal_node_ids:
        node_roles[node_i] = 2

    graph_edges: list[tuple[int, int]] = []
    edge_original_ids: list[int] = []

    for edge_i, (start_i, end_i) in enumerate(edges):
        start = nodes[start_i]
        end = nodes[end_i]
        length = math.dist(start, end)
        n_segments = max(1, math.ceil(length / graph_step))

        previous_i = start_i
        for segment_i in range(1, n_segments):
            alpha = segment_i / n_segments
            point = tuple(
                (1.0 - alpha) * start[component_i]
                + alpha * end[component_i]
                for component_i in range(3)
            )
            point_i = len(graph_points)
            graph_points.append(point)
            node_roles.append(0)
            graph_edges.append((previous_i, point_i))
            edge_original_ids.append(edge_i)
            previous_i = point_i

        graph_edges.append((previous_i, end_i))
        edge_original_ids.append(edge_i)

    print(f"Writing line graph VTK: {graph_path}")
    with graph_path.open("w", encoding="utf-8") as vtk:
        vtk.write("# vtk DataFile Version 2.0\n")
        vtk.write("Purkinje graph for 1DgraphToFoam\n")
        vtk.write("ASCII\n")
        vtk.write("DATASET UNSTRUCTURED_GRID\n")
        vtk.write(f"POINTS {len(graph_points)} double\n")
        for point in graph_points:
            vtk.write(f"{point[0]:.12g} {point[1]:.12g} {point[2]:.12g}\n")

        vtk.write(f"CELLS {len(graph_edges)} {3 * len(graph_edges)}\n")
        for start_i, end_i in graph_edges:
            vtk.write(f"2 {start_i} {end_i}\n")

        vtk.write(f"CELL_TYPES {len(graph_edges)}\n")
        for _ in graph_edges:
            vtk.write("3\n")  # VTK_LINE

        vtk.write(f"POINT_DATA {len(graph_points)}\n")
        vtk.write("SCALARS nodeRole int 1\n")
        vtk.write("LOOKUP_TABLE default\n")
        for role in node_roles:
            vtk.write(f"{role}\n")

        vtk.write(f"CELL_DATA {len(graph_edges)}\n")
        vtk.write("SCALARS originalEdgeId int 1\n")
        vtk.write("LOOKUP_TABLE default\n")
        for edge_i in edge_original_ids:
            vtk.write(f"{edge_i}\n")

        vtk.write("SCALARS radius double 1\n")
        vtk.write("LOOKUP_TABLE default\n")
        for _ in graph_edges:
            vtk.write(f"{radius:.12g}\n")

        vtk.write("SCALARS diffusivity double 1\n")
        vtk.write("LOOKUP_TABLE default\n")
        for _ in graph_edges:
            vtk.write(f"{diffusivity:.12g}\n")

        vtk.write("SCALARS conductance double 1\n")
        vtk.write("LOOKUP_TABLE default\n")
        for _ in graph_edges:
            vtk.write(f"{diffusivity:.12g}\n")


def write_outputs(
    output_dir: Path,
    basename: str,
    graph_step: float,
    radius: float,
    diffusivity: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    msh_path = output_dir / f"{basename}.msh"
    vtk_path = output_dir / f"{basename}.vtk"

    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)

    # gmshToFoam is most predictable with legacy ASCII MSH 2.2.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    print(f"Writing OpenFOAM conversion mesh: {msh_path}")
    gmsh.write(str(msh_path))

    # Diagnostic file: useful to inspect CellEntityIds/physical groups in ParaView.
    gmsh.option.setNumber("Mesh.Format", 42)
    print(f"Writing diagnostic VTK: {vtk_path}")
    gmsh.write(str(vtk_path))

    write_graph_vtk(
        output_dir,
        basename,
        NODES,
        EDGES,
        TERMINAL_NODE_IDS,
        radius,
        graph_step,
        diffusivity,
    )


def main() -> None:
    args = parse_args()

    if args.graph_only:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_graph_vtk(
            args.output_dir,
            args.basename,
            NODES,
            EDGES,
            TERMINAL_NODE_IDS,
            args.radius,
            args.graph_step,
            args.diffusivity,
        )
        print("Done. Use 1DgraphToFoam on the *_graph.vtk file.")
        return

    global gmsh
    import gmsh

    gmsh.initialize()
    try:
        gmsh.model.add("PurkinjeNetwork")

        add_tube_volume(NODES, EDGES, args.radius)
        add_physical_groups(NODES, args.radius)

        gmsh.option.setNumber("Mesh.MeshSizeMin", args.mesh_size * 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", args.mesh_size)

        print("Generating 3D Purkinje mesh...")
        gmsh.model.mesh.generate(3)

        write_outputs(
            args.output_dir,
            args.basename,
            args.graph_step,
            args.radius,
            args.diffusivity,
        )
    finally:
        gmsh.finalize()

    print("Done. Use 1DgraphToFoam on the *_graph.vtk file.")


if __name__ == "__main__":
    main()
