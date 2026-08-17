#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     plugins.cardiacfoam.mesh_geometry
#
# Description
#     Plan-time SI-scale checks for purkinjeGraph dictionaries, layered on
#     core's generic polyMesh bounding-box machinery.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""purkinjeGraph scale detection (cardiacFoam's point sets that are not meshes).

`specs/mesh_geometry.py` classifies the scale of every polyMesh region in a
case. A cardiacFoam case can also carry a Purkinje conduction tree in
`constant/purkinjeGraph*`, which is a Foam *dictionary* holding its own point
list -- not a mesh region, so core's region discovery never sees it, yet it
must share the mesh's coordinate frame or the PVJ coupling lands nowhere.

This module applies core's own `classify_scale` to those graphs and
cross-checks them against the default mesh region, and is wired into the
strict planner as the cardiac plugin's extra mesh-geometry diagnostics.
Stdlib-only; never mutates the case.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from openfoam_driver.specs.mesh_geometry import (
    ASCII_TRIPLE_RE,
    BoundingBox,
    MeshDiagnostic,
    MeshParseError,
    ScaleClass,
    bounding_box_from_flat_coords,
    classify_scale,
    discover_mesh_regions,
    read_bounding_box,
)

# Matches the top-level `points` section in a Foam dictionary (e.g. purkinjeGraph).
# \b prevents matching `pointFields` or `endpointNodes`.
_GRAPH_POINTS_RE = re.compile(rb"\bpoints\b\s+(\d+)\s*\(")


def read_purkinje_graph_bbox(graph_path: Path) -> BoundingBox:
    """Return the bounding box of the points section in a purkinjeGraph dict.

    purkinjeGraph is a Foam dictionary (not a standalone field file), so we
    locate the ``points`` keyword explicitly rather than reusing
    read_bounding_box, which would land on the earlier ``pvjNodes`` integer
    list and find no coordinate triples.
    """
    data = graph_path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)

    m = _GRAPH_POINTS_RE.search(data)
    if not m:
        raise MeshParseError(
            f"could not locate 'points' section in {graph_path.name}"
        )
    count = int(m.group(1))
    if count == 0:
        raise MeshParseError(f"{graph_path.name} declares zero graph points")

    text = data[m.end():].decode("latin-1")
    triples = ASCII_TRIPLE_RE.findall(text)
    if not triples:
        raise MeshParseError(f"no coordinate triples found after 'points' in {graph_path.name}")

    coords = [float(v) for triple in triples for v in triple]
    return bounding_box_from_flat_coords(coords)


def discover_purkinje_graphs(case_root: Path) -> list[Path]:
    """Find all purkinjeGraph* dictionary files under constant/.

    Returns plain files (not directories) whose names start with
    ``purkinjeGraph``.  Typical examples: ``purkinjeGraph``,
    ``purkinjeGraphScar``.
    """
    constant = case_root / "constant"
    if not constant.is_dir():
        return []
    return sorted(
        p for p in constant.iterdir()
        if p.is_file() and p.name.startswith("purkinjeGraph")
    )


def _default_region_scale(case_root: Path) -> ScaleClass | None:
    """Scale of the default (unnamed) mesh region, or None if unavailable.

    A graph is only cross-checked against the mesh when the default region
    parsed cleanly -- an unparseable or missing points file already produces
    its own core diagnostic, and guessing a unit from it would turn one
    problem into two.
    """
    for region in discover_mesh_regions(case_root):
        if region.name:
            continue
        try:
            return classify_scale(read_bounding_box(region.points_path).max_dim)
        except (MeshParseError, OSError):
            return None
    return None


def purkinje_graph_diagnostics(case_root: Path) -> tuple[MeshDiagnostic, ...]:
    """Flag non-SI Purkinje graphs and graph/mesh scale disagreement.

    Same `classify_scale` logic core applies to mesh regions, cross-checked
    against the default mesh region when available.

    Like core's own mesh-scale gate, this returns nothing for a case with no
    mesh region at all: the graph check has always been a companion to the
    mesh check, and a case with no mesh has bigger problems that the planner
    reports elsewhere. Preserved deliberately from when this loop lived inside
    `specs/mesh_geometry.py::mesh_geometry_diagnostics`, whose
    ``if not regions: return ()`` guard short-circuited it.
    """
    if not discover_mesh_regions(case_root):
        return ()

    graphs = discover_purkinje_graphs(case_root)
    if not graphs:
        return ()

    mesh_scale = _default_region_scale(case_root)
    diagnostics: list[MeshDiagnostic] = []

    for graph_path in graphs:
        name = graph_path.name
        try:
            gbbox = read_purkinje_graph_bbox(graph_path)
        except (MeshParseError, OSError) as exc:
            diagnostics.append(MeshDiagnostic(
                "warning", "graph_scale_not_checked",
                f"Could not parse points from '{name}': {exc}. "
                f"Verify that the graph was generated in SI metres.",
                name,
            ))
            continue
        gscale = classify_scale(gbbox.max_dim)
        if gscale.unit != "m":
            diagnostics.append(MeshDiagnostic(
                "error", "graph_not_si",
                f"'{name}' max dimension {gbbox.max_dim:g} suggests "
                f"{gscale.unit}, not metres. Regenerate the Purkinje tree "
                f"after rescaling the mesh with 'checkMeshGeometry -rescale'.",
                name,
            ))
        if mesh_scale is not None:
            if gscale.unit != mesh_scale.unit:
                diagnostics.append(MeshDiagnostic(
                    "error", "graph_mesh_scale_mismatch",
                    f"'{name}' is in {gscale.unit} but the mesh is in "
                    f"{mesh_scale.unit}. Rescale both to SI metres before running.",
                    name,
                ))

    return tuple(diagnostics)
