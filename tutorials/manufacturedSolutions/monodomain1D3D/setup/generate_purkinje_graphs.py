#!/usr/bin/env python3
"""Generate refined Y-shaped Purkinje graph dictionaries and VTK views."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


DEFAULT_GRAPH_Y = 1.0 / 6.0
DEFAULT_GRAPH_Z = 1.0 / 3.0
CONDUCTANCE = 1.0


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def interpolate(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    fraction: float,
) -> tuple[float, float, float]:
    return tuple(ai + fraction * (bi - ai) for ai, bi in zip(a, b))


def fmt_point(point: tuple[float, float, float]) -> str:
    return f"({point[0]:.12g} {point[1]:.12g} {point[2]:.12g})"


def build_graph(segments_per_branch: int, graph_y: float, graph_z: float):
    if segments_per_branch < 1:
        raise ValueError("segments_per_branch must be positive")

    root = (0.50, graph_y, graph_z)
    terminals = (
        (0.00, graph_y, graph_z),
        (1.00, graph_y, graph_z),
    )

    points = [root]
    edges = []
    pvj_nodes = []

    for terminal in terminals:
        previous = 0
        edge_length = distance(root, terminal) / segments_per_branch

        for segment_i in range(1, segments_per_branch + 1):
            points.append(interpolate(root, terminal, segment_i / segments_per_branch))
            node = len(points) - 1
            edges.append((previous, node, edge_length, CONDUCTANCE))
            previous = node

        pvj_nodes.append(previous)

    return points, edges, pvj_nodes


def write_foam_graph(path: Path, points, edges, pvj_nodes) -> None:
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      purkinjeGraph;
}

rootNode
0;

pvjNodes
"""
        )
        stream.write("(" + " ".join(str(node) for node in pvj_nodes) + ");\n\n")
        stream.write("pvjLocations\n(\n")
        for node in pvj_nodes:
            stream.write(f"    {fmt_point(points[node])}\n")
        stream.write(");\n\n")

        stream.write("conductionEdges\n(\n")
        for a, b, length, conductance in edges:
            stream.write(f"    ({a} {b} {length:.12g} {conductance:.12g})\n")
        stream.write(");\n\n")

        stream.write("points\n(\n")
        for point in points:
            stream.write(f"    {fmt_point(point)}\n")
        stream.write(");\n")


def write_vtk_graph(path: Path, points, edges, pvj_nodes) -> None:
    pvj_set = set(pvj_nodes)

    with path.open("w", encoding="utf-8") as stream:
        stream.write("# vtk DataFile Version 3.0\n")
        stream.write(f"Purkinje graph with {len(points)} nodes\n")
        stream.write("ASCII\n")
        stream.write("DATASET POLYDATA\n")

        stream.write(f"POINTS {len(points)} float\n")
        for point in points:
            stream.write(f"{point[0]:.12g} {point[1]:.12g} {point[2]:.12g}\n")

        stream.write(f"LINES {len(edges)} {len(edges) * 3}\n")
        for a, b, _, _ in edges:
            stream.write(f"2 {a} {b}\n")

        stream.write(f"POINT_DATA {len(points)}\n")
        stream.write("SCALARS nodeRole int 1\n")
        stream.write("LOOKUP_TABLE default\n")
        for node_i in range(len(points)):
            if node_i == 0:
                role = 1
            elif node_i in pvj_set:
                role = 2
            else:
                role = 0
            stream.write(f"{role}\n")

        stream.write(f"CELL_DATA {len(edges)}\n")
        stream.write("SCALARS edgeLength float 1\n")
        stream.write("LOOKUP_TABLE default\n")
        for _, _, length, _ in edges:
            stream.write(f"{length:.12g}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Tutorial case directory",
    )
    parser.add_argument(
        "--segments",
        type=int,
        nargs="+",
        default=(1, 5, 10, 20, 40, 80),
        help="Segment counts per branch",
    )
    parser.add_argument(
        "--graph-y",
        type=float,
        default=DEFAULT_GRAPH_Y,
        help="Constant y-coordinate for the 1D graph",
    )
    parser.add_argument(
        "--graph-z",
        type=float,
        default=DEFAULT_GRAPH_Z,
        help="Constant z-coordinate for the 1D graph",
    )
    args = parser.parse_args()

    case_dir = args.case.resolve()
    constant_dir = case_dir / "constant"
    vtk_dir = constant_dir / "graphFiles"
    vtk_dir.mkdir(parents=True, exist_ok=True)

    for segments in args.segments:
        points, edges, pvj_nodes = build_graph(segments, args.graph_y, args.graph_z)
        node_count = len(points)
        suffix = f"nodes{node_count:03d}"
        foam_path = constant_dir / f"purkinjeGraph.{suffix}"
        vtk_path = vtk_dir / f"purkinjeGraph.{suffix}.vtk"

        write_foam_graph(foam_path, points, edges, pvj_nodes)
        write_vtk_graph(vtk_path, points, edges, pvj_nodes)

        print(
            f"{suffix}: segments/branch={segments}, "
            f"nodes={len(points)}, edges={len(edges)}, "
            f"graph_y={args.graph_y:.12g}, graph_z={args.graph_z:.12g}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
