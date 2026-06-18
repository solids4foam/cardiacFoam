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
#     test_mesh_geometry
#
# Description
#     Unit tests for the plan-time mesh-scale detection module.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import gzip
import struct
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.specs.mesh_geometry import (
    BoundingBox,
    MeshDiagnostic,
    MeshParseError,
    MeshRegion,
    ScaleClass,
    classify_scale,
    discover_mesh_regions,
    mesh_geometry_diagnostics,
    read_bounding_box,
)
from openfoam_driver.utility_catalog import UTILITY_CATALOG


_HEADER = (
    "FoamFile\n{{\n    version 2.0;\n    format {fmt};\n"
    '    arch "LSB;label=32;scalar=64";\n'
    "    class vectorField;\n    object points;\n}}\n"
    "// * * * //\n"
)


def _write_ascii_points(path: Path, pts):
    body = "".join(f"({x} {y} {z})\n" for x, y, z in pts)
    path.write_text(_HEADER.format(fmt="ascii") + f"\n{len(pts)}\n(\n{body})\n")


def _write_binary_points(path: Path, pts):
    # Real OpenFOAM emits a newline after the opening '(' before the binary
    # block; the parser must skip it. Reproduce that here.
    header = _HEADER.format(fmt="binary").encode("latin-1")
    flat = [c for triple in pts for c in triple]
    block = struct.pack(f"<{len(flat)}d", *flat)
    path.write_bytes(header + f"\n{len(pts)}\n(\n".encode("latin-1") + block + b"\n)\n")


class TestClassifyScale(unittest.TestCase):
    def test_metres_passthrough(self):
        self.assertEqual(classify_scale(0.12), ScaleClass("m", 1.0))

    def test_millimetre_range(self):
        self.assertEqual(classify_scale(50.0), ScaleClass("mm", 1e-3))

    def test_millimetre_cutoff_inclusive(self):
        # 20.0 is the metres/mm cutoff (inclusive on the mm side).
        self.assertEqual(classify_scale(20.0), ScaleClass("mm", 1e-3))

    def test_large_si_domain_is_metres(self):
        # A whole-torso / whole-body SI mesh (~1.7 m) must NOT be flagged mm.
        self.assertEqual(classify_scale(1.7), ScaleClass("m", 1.0))
        self.assertEqual(classify_scale(19.9), ScaleClass("m", 1.0))

    def test_micrometre_range(self):
        self.assertEqual(classify_scale(5000.0), ScaleClass("um", 1e-6))

    def test_micrometre_lower_bound_inclusive(self):
        self.assertEqual(classify_scale(1000.0), ScaleClass("um", 1e-6))

    def test_very_large_falls_back_to_metres(self):
        self.assertEqual(classify_scale(2e6), ScaleClass("m", 1.0))


class TestReadBoundingBox(unittest.TestCase):
    PTS = [(0.0, 0.0, 0.0), (2.0, 5.0, 1.0), (-1.0, 3.0, 4.0)]

    def test_ascii(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "points"
            _write_ascii_points(p, self.PTS)
            bb = read_bounding_box(p)
            self.assertEqual(bb, BoundingBox((-1.0, 0.0, 0.0), (2.0, 5.0, 4.0)))
            self.assertEqual(bb.max_dim, 5.0)

    def test_binary(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "points"
            _write_binary_points(p, self.PTS)
            bb = read_bounding_box(p)
            self.assertEqual(bb, BoundingBox((-1.0, 0.0, 0.0), (2.0, 5.0, 4.0)))

    def test_gzip_binary(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d) / "points"
            _write_binary_points(raw, self.PTS)
            gz = Path(d) / "points.gz"
            gz.write_bytes(gzip.compress(raw.read_bytes()))
            bb = read_bounding_box(gz)
            self.assertEqual(bb.max_dim, 5.0)

    def test_truncated_binary_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "points"
            header = _HEADER.format(fmt="binary").encode("latin-1")
            p.write_bytes(header + b"\n3\n(" + struct.pack("<3d", 1.0, 2.0, 3.0))
            with self.assertRaises(MeshParseError):
                read_bounding_box(p)


class TestDiscoverMeshRegions(unittest.TestCase):
    def test_default_region(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            pm = case / "constant" / "polyMesh"
            pm.mkdir(parents=True)
            _write_ascii_points(pm / "points", TestReadBoundingBox.PTS)
            regions = discover_mesh_regions(case)
            self.assertEqual([r.name for r in regions], [""])

    def test_multi_region_sorted(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            for name in ("myocardium", "purkinje"):
                pm = case / "constant" / name / "polyMesh"
                pm.mkdir(parents=True)
                _write_ascii_points(pm / "points", TestReadBoundingBox.PTS)
            regions = discover_mesh_regions(case)
            self.assertEqual([r.name for r in regions], ["myocardium", "purkinje"])

    def test_gz_points_discovered(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            pm = case / "constant" / "polyMesh"
            pm.mkdir(parents=True)
            raw = pm / "points"
            _write_ascii_points(raw, TestReadBoundingBox.PTS)
            (pm / "points.gz").write_bytes(gzip.compress(raw.read_bytes()))
            raw.unlink()
            regions = discover_mesh_regions(case)
            self.assertEqual([r.name for r in regions], [""])
            self.assertTrue(regions[0].points_path.name.endswith(".gz"))

    def test_no_mesh_returns_empty(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(discover_mesh_regions(Path(d)), [])


def _make_case(d, region_pts: dict, fmt="ascii"):
    """region_pts: {region_name_or_'': [(x,y,z), ...]}."""
    case = Path(d)
    writer = _write_ascii_points if fmt == "ascii" else _write_binary_points
    for name, pts in region_pts.items():
        sub = "polyMesh" if name == "" else f"{name}/polyMesh"
        pm = case / "constant" / sub
        pm.mkdir(parents=True)
        writer(pm / "points", pts)
    return case


_SI = [(0.0, 0.0, 0.0), (0.05, 0.05, 0.05)]       # max_dim 0.05 m
_MM = [(0.0, 0.0, 0.0), (50.0, 50.0, 50.0)]        # max_dim 50 (>=20) -> mm
_FAR = [(100.0, 0.0, 0.0), (100.05, 0.05, 0.05)]   # SI but disjoint bbox
_TORSO = [(0.0, 0.0, 0.0), (1.7, 1.7, 1.7)]        # max_dim 1.7 -> m, ambiguous


class TestMeshGeometryDiagnostics(unittest.TestCase):
    def test_si_mesh_is_clean(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case(d, {"": _SI})
            self.assertEqual(mesh_geometry_diagnostics(case), ())

    def test_non_si_mesh_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case(d, {"": _MM})
            codes = {x.code for x in mesh_geometry_diagnostics(case)}
            self.assertIn("mesh_not_si", codes)

    def test_large_si_domain_is_ambiguous_warning_not_error(self):
        # A whole-torso SI mesh (~1.7 m) must NOT be flagged mm; it gets an
        # advisory warning instead, so the agent never tries to "rescale" it.
        with tempfile.TemporaryDirectory() as d:
            case = _make_case(d, {"": _TORSO})
            diags = mesh_geometry_diagnostics(case)
            codes = {x.code for x in diags}
            self.assertNotIn("mesh_not_si", codes)
            self.assertIn("mesh_scale_ambiguous", codes)
            warn = [x for x in diags if x.code == "mesh_scale_ambiguous"][0]
            self.assertEqual(warn.level, "warning")

    def test_coupled_regions_unit_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case(d, {"myocardium": _SI, "purkinje": _MM})
            diags = mesh_geometry_diagnostics(case)
            codes = {x.code for x in diags}
            self.assertIn("mesh_not_si", codes)          # purkinje is mm
            self.assertIn("mesh_scale_mismatch", codes)  # regions disagree

    def test_same_unit_disjoint_bbox_mismatch_is_warning(self):
        with tempfile.TemporaryDirectory() as d:
            case = _make_case(d, {"myocardium": _SI, "purkinje": _FAR})
            diags = mesh_geometry_diagnostics(case)
            mismatch = [x for x in diags if x.code == "mesh_scale_mismatch"]
            self.assertTrue(mismatch)
            # same-unit non-overlap is advisory, not a hard block
            self.assertEqual(mismatch[0].level, "warning")

    def test_unparseable_points_is_loud_not_silent(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            pm = case / "constant" / "polyMesh"
            pm.mkdir(parents=True)
            (pm / "points").write_text("not a points file")
            diags = mesh_geometry_diagnostics(case)
            self.assertEqual([x.code for x in diags], ["mesh_scale_not_checked"])
            self.assertEqual(diags[0].level, "warning")

    def test_no_mesh_no_diagnostics(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(mesh_geometry_diagnostics(Path(d)), ())


class TestCheckMeshGeometryCatalogued(unittest.TestCase):
    def test_present_and_mesh_category(self):
        self.assertIn("checkMeshGeometry", UTILITY_CATALOG)
        entry = UTILITY_CATALOG["checkMeshGeometry"]
        self.assertEqual(entry.category, "mesh")
        self.assertTrue(entry.requires_mesh)

    def test_flags_reflect_detect_only_default(self):
        flag_names = {f.name for f in UTILITY_CATALOG["checkMeshGeometry"].flags}
        self.assertIn("-region", flag_names)
        self.assertIn("-scale", flag_names)
        self.assertIn("-rescale", flag_names)
        # -noScale is gone: detect-only is now the default.
        self.assertNotIn("-noScale", flag_names)

    def test_artifact_id_typo_fixed(self):
        produced = {p.artifact_id for p in UTILITY_CATALOG["checkMeshGeometry"].produces}
        self.assertIn("polymesh_scaled", produced)
        self.assertNotIn("polymes_scaled", produced)


if __name__ == "__main__":
    unittest.main()
