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
#     test_purkinje_graph_geometry
#
# Description
#     Tests the cardiac plugin's purkinjeGraph SI-scale diagnostics, layered
#     on core's generic polyMesh bounding-box machinery.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.plugins.cardiacfoam.mesh_geometry import (
    discover_purkinje_graphs,
    purkinje_graph_diagnostics,
    read_purkinje_graph_bbox,
)
from openfoam_driver.core.specs.mesh_geometry import MeshParseError


_HEADER = (
    "FoamFile\n{{\n    version 2.0;\n    format {fmt};\n"
    '    arch "LSB;label=32;scalar=64";\n'
    "    class vectorField;\n    object points;\n}}\n"
    "// * * * //\n"
)

_GRAPH_HEADER = (
    "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
    "    class dictionary;\n    object purkinjeGraph;\n}\n\n"
)


def _write_ascii_points(path: Path, pts):
    body = "".join(f"({x} {y} {z})\n" for x, y, z in pts)
    path.write_text(_HEADER.format(fmt="ascii") + f"\n{len(pts)}\n(\n{body})\n")


def _write_purkinje_graph(path: Path, pts, pvj_pts=None):
    """Write a minimal purkinjeGraph Foam dict with the given points section."""
    if pvj_pts is None:
        pvj_pts = pts[:1]
    n_pvj = len(pvj_pts)
    n_pts = len(pts)
    # pvjNodes (integer indices, no triples)
    pvj_nodes = "".join(f"{i}\n" for i in range(n_pvj))
    pvj_locs = "".join(f"({x} {y} {z})\n" for x, y, z in pvj_pts)
    pt_body = "".join(f"({x} {y} {z})\n" for x, y, z in pts)
    content = (
        _GRAPH_HEADER
        + f"rootNode\n0;\n\n"
        + f"pvjNodes\n\n{n_pvj}\n(\n{pvj_nodes})\n;\n\n"
        + f"pvjLocations\n\n{n_pvj}\n(\n{pvj_locs})\n;\n\n"
        + f"conductionEdges\n\n0\n(\n)\n;\n\n"
        + f"points\n\n{n_pts}\n(\n{pt_body})\n;\n\n"
    )
    path.write_text(content)


def _make_case_with_graph(d, mesh_pts, graph_pts, graph_name="purkinjeGraph"):
    case = Path(d)
    pm = case / "constant" / "polyMesh"
    pm.mkdir(parents=True)
    _write_ascii_points(pm / "points", mesh_pts)
    _write_purkinje_graph(case / "constant" / graph_name, graph_pts)
    return case


_SI = [(0.0, 0.0, 0.0), (0.05, 0.05, 0.05)]       # max_dim 0.05 m
_MM = [(0.0, 0.0, 0.0), (50.0, 50.0, 50.0)]        # max_dim 50 (>=20) -> mm
_MM_GRAPH = [(-5.0, 80.0, -200.0), (40.0, 110.0, -170.0)]   # max_dim ~45 -> mm
_SI_GRAPH = [(0.0, 0.0, 0.0), (0.04, 0.11, 0.07)]            # max_dim 0.11 -> m


class TestReadPurkinjeGraphBbox(unittest.TestCase):
    PTS = [(-5.0, 80.0, -200.0), (40.0, 110.0, -170.0)]  # mm-scale biventricular

    def test_reads_points_section(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "purkinjeGraph"
            _write_purkinje_graph(p, self.PTS)
            bb = read_purkinje_graph_bbox(p)
            self.assertAlmostEqual(bb.max_dim, 45.0)  # X span: -5 to 40

    def test_pvj_nodes_integers_do_not_corrupt_bbox(self):
        # pvjNodes contains bare integers (not triples); they must not be
        # mistaken for coordinates. The bbox must reflect only the points section.
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "purkinjeGraph"
            _write_purkinje_graph(p, self.PTS, pvj_pts=[(10.0, 90.0, -190.0)])
            bb = read_purkinje_graph_bbox(p)
            # X range of points only: -5 to 40 → span 45
            self.assertAlmostEqual(bb.min_pt[0], -5.0)
            self.assertAlmostEqual(bb.max_pt[0], 40.0)

    def test_missing_points_section_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "purkinjeGraph"
            p.write_text(_GRAPH_HEADER + "rootNode\n0;\n")
            with self.assertRaises(MeshParseError):
                read_purkinje_graph_bbox(p)

    def test_si_scale_points(self):
        si_pts = [(0.0, 0.0, 0.0), (0.04, 0.11, 0.07)]  # metres
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "purkinjeGraph"
            _write_purkinje_graph(p, si_pts)
            bb = read_purkinje_graph_bbox(p)
            self.assertLess(bb.max_dim, 20.0)


class TestDiscoverPurkinjeGraphs(unittest.TestCase):
    def test_finds_purkinjeGraph(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            (case / "constant").mkdir()
            p = case / "constant" / "purkinjeGraph"
            _write_purkinje_graph(p, [(0, 0, 0)])
            found = discover_purkinje_graphs(case)
            self.assertEqual([f.name for f in found], ["purkinjeGraph"])

    def test_finds_scar_variant(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            (case / "constant").mkdir()
            for name in ("purkinjeGraph", "purkinjeGraphScar"):
                _write_purkinje_graph(case / "constant" / name, [(0, 0, 0)])
            found = discover_purkinje_graphs(case)
            self.assertEqual([f.name for f in found], ["purkinjeGraph", "purkinjeGraphScar"])

    def test_ignores_directories_named_purkinjeGraph(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            (case / "constant" / "purkinjeGraph").mkdir(parents=True)
            found = discover_purkinje_graphs(case)
            self.assertEqual(found, [])

    def test_no_constant_dir_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(discover_purkinje_graphs(Path(d)), [])


class TestPurkinjeGraphDiagnostics(unittest.TestCase):
    def test_si_graph_with_si_mesh_is_clean(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _SI, _SI_GRAPH)
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertNotIn("graph_not_si", codes)
            self.assertNotIn("graph_mesh_scale_mismatch", codes)

    def test_mm_graph_flagged_as_not_si(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _SI, _MM_GRAPH)
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertIn("graph_not_si", codes)

    def test_mm_graph_with_mm_mesh_flagged_graph_only_no_mismatch(self):
        # Both mesh and graph are in mm: the graph is flagged non-SI, but the
        # two agree, so no mismatch. (mesh_not_si is core's own diagnostic.)
        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _MM, _MM_GRAPH)
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertIn("graph_not_si", codes)
            self.assertNotIn("graph_mesh_scale_mismatch", codes)

    def test_scale_mismatch_between_graph_and_mesh(self):
        # SI mesh but mm graph → mismatch error in addition to graph_not_si.
        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _SI, _MM_GRAPH)
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertIn("graph_not_si", codes)
            self.assertIn("graph_mesh_scale_mismatch", codes)

    def test_scar_graph_also_checked(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _SI, _MM_GRAPH, graph_name="purkinjeGraphScar")
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertIn("graph_not_si", codes)

    def test_unparseable_graph_is_warning_not_error(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            pm = case / "constant" / "polyMesh"
            pm.mkdir(parents=True)
            _write_ascii_points(pm / "points", _SI)
            (case / "constant" / "purkinjeGraph").write_text("not a graph file")
            codes_levels = {
                (x.code, x.level) for x in purkinje_graph_diagnostics(case)
            }
            self.assertIn(("graph_scale_not_checked", "warning"), codes_levels)

    def test_unparseable_mesh_suppresses_only_the_cross_check(self):
        # The mesh's own unit is unknown, so no mismatch can be asserted --
        # but the graph is still classified on its own.
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            pm = case / "constant" / "polyMesh"
            pm.mkdir(parents=True)
            (pm / "points").write_text("not a points file")
            _write_purkinje_graph(case / "constant" / "purkinjeGraph", _MM_GRAPH)
            codes = {x.code for x in purkinje_graph_diagnostics(case)}
            self.assertIn("graph_not_si", codes)
            self.assertNotIn("graph_mesh_scale_mismatch", codes)

    def test_no_mesh_region_means_no_graph_diagnostics(self):
        # Preserved from when this loop lived inside mesh_geometry_diagnostics,
        # whose `if not regions: return ()` guard short-circuited it.
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            (case / "constant").mkdir()
            _write_purkinje_graph(case / "constant" / "purkinjeGraph", _MM_GRAPH)
            self.assertEqual(purkinje_graph_diagnostics(case), ())


class TestPurkinjeDiagnosticsReachTheStrictPlanner(unittest.TestCase):
    def test_plugin_hook_is_wired_into_mesh_geometry_diagnostics(self):
        # The move only holds if the plugin's checks still reach the report.
        from openfoam_driver.core.strict_planning import _mesh_geometry_diagnostics

        with tempfile.TemporaryDirectory() as d:
            case = _make_case_with_graph(d, _SI, _MM_GRAPH)
            diags = _mesh_geometry_diagnostics(case)
            codes = {x.code for x in diags}
            self.assertIn("graph_not_si", codes)
            self.assertTrue(all(x.source == "mesh_geometry" for x in diags))


if __name__ == "__main__":
    unittest.main()
