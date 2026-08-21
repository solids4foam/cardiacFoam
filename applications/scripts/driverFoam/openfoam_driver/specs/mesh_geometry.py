#----------------------------------------------------------------------------#
# Module
#     mesh_geometry
#
# Description
#     Non-mutating, plan-time detection of mesh point scale. Reads polyMesh
#     point bounding boxes directly (ASCII/binary/gz) so the strict planner can
#     flag non-SI meshes and coupled-region scale mismatches before any compute.
#     Plugins can add their own point-set checks (e.g. a conduction graph) on
#     top, reusing this module's parsing and classification primitives.
#     The rescaling *write* is delegated to the checkMeshGeometry utility.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#
"""Plan-time mesh-scale detection. Stdlib-only; never mutates the case."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from foamlib import FoamFile


# Thresholds MUST mirror applications/utilities/checkMeshGeometry/checkMeshGeometry.C.
# A drift guard lives in test_mesh_geometry_contract.py.
#
# _MM_LOWER is the metres/mm cutoff. It is 20.0 (not 1.0) so that large SI
# domains -- whole-torso / whole-body bidomain meshes (~1-2 m) -- are NOT
# mis-flagged as millimetres, while a heart authored in mm (~50-150) still
# lands at/above the cutoff and is flagged. This is the one knob to tune if
# your physical meshes are larger than ~20 m in SI (they should not be).
_MM_LOWER = 20.0
_UM_LOWER = 1000.0
_UM_UPPER = 1.0e6

# Below the mm cutoff but >= 1.0 raw units is ambiguous: it could be a large SI
# domain (a whole-torso mesh) OR a small mm/cm mesh. The gate treats this band
# as metres (no rescale) but emits an advisory warning -- it never hard-blocks.
_AMBIGUOUS_LOWER = 1.0


@dataclass(frozen=True)
class ScaleClass:
    """Detected unit and the factor that converts the mesh to SI metres."""

    unit: str
    scale_factor: float


def classify_scale(max_dim: float) -> ScaleClass:
    """Classify a mesh by its largest bounding-box extent (raw units).

    max_dim < 20      -> "m"  (covers all SI cardiac/torso/whole-body domains)
    20  <= max_dim < 1e3 -> "mm"
    1e3 <= max_dim < 1e6 -> "um"
    max_dim >= 1e6    -> "m"  (implausibly large; fall back to metres)
    """
    if _MM_LOWER <= max_dim < _UM_LOWER:
        return ScaleClass("mm", 1e-3)
    if _UM_LOWER <= max_dim < _UM_UPPER:
        return ScaleClass("um", 1e-6)
    return ScaleClass("m", 1.0)


class MeshParseError(Exception):
    """Raised when a points file exists but cannot be parsed for a bbox."""


@dataclass(frozen=True)
class BoundingBox:
    min_pt: tuple[float, float, float]
    max_pt: tuple[float, float, float]

    @property
    def max_dim(self) -> float:
        return max(hi - lo for lo, hi in zip(self.min_pt, self.max_pt))


# One `(x y z)` vector in an OpenFOAM ASCII list. Public so a plugin adding
# its own point-set check parses coordinates exactly the way core does.
ASCII_TRIPLE_RE = re.compile(
    r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)"
)


def bounding_box_from_flat_coords(coords) -> BoundingBox:
    """Bounding box of a flat ``[x0, y0, z0, x1, ...]`` coordinate sequence."""
    if len(coords) < 3:
        raise MeshParseError("fewer than one point parsed")
    xs, ys, zs = coords[0::3], coords[1::3], coords[2::3]
    return BoundingBox(
        (min(xs), min(ys), min(zs)),
        (max(xs), max(ys), max(zs)),
    )


def read_bounding_box(points_path: Path) -> BoundingBox:
    """Return the point bounding box, handling ASCII/binary/gz formats."""
    try:
        points = FoamFile(points_path)[None]
        if len(points) == 0:
            raise MeshParseError("fewer than one point parsed")
        coords = [float(v) for point in points for v in point]
    except (ValueError, OSError, TypeError) as exc:
        # Not just a foamlib decode failure: a syntactically valid but
        # non-points file (e.g. bare standalone words) parses without
        # raising, so the failure only surfaces once we try to read the
        # result as coordinates. Both cases map to the same MeshParseError.
        raise MeshParseError(f"could not parse points file: {exc}") from exc
    return bounding_box_from_flat_coords(coords)


@dataclass(frozen=True)
class MeshRegion:
    name: str  # "" for the default region
    points_path: Path


def _points_file(polymesh_dir: Path) -> Path | None:
    for candidate in ("points", "points.gz"):
        path = polymesh_dir / candidate
        if path.exists():
            return path
    return None


def discover_mesh_regions(case_root: Path) -> list[MeshRegion]:
    """Find every region with a polyMesh under constant/ (default + named)."""
    regions: list[MeshRegion] = []
    constant = case_root / "constant"
    if not constant.is_dir():
        return regions

    default_pts = _points_file(constant / "polyMesh")
    if default_pts is not None:
        regions.append(MeshRegion("", default_pts))

    for child in sorted(constant.iterdir()):
        if not child.is_dir() or child.name == "polyMesh":
            continue
        pts = _points_file(child / "polyMesh")
        if pts is not None:
            regions.append(MeshRegion(child.name, pts))

    return regions


@dataclass(frozen=True)
class MeshDiagnostic:
    level: str   # "error" | "warning"
    code: str
    message: str
    region: str = ""


def _bbox_overlaps(a: BoundingBox, b: BoundingBox) -> bool:
    return all(
        a.min_pt[i] <= b.max_pt[i] and b.min_pt[i] <= a.max_pt[i]
        for i in range(3)
    )


def mesh_geometry_diagnostics(
    case_root: Path,
    *,
    coupled_groups: list[frozenset[str]] | None = None,
) -> tuple[MeshDiagnostic, ...]:
    """Detect non-SI meshes and coupled-region scale disagreement.

    coupled_groups names regions that must share scale. When omitted, all
    discovered regions are treated as one implicitly-coupled group — a
    conservative heuristic that may over-warn for regions that are not actually
    coupled. Deriving real groups from the case's own declared region couplings
    is a deferred follow-up. To limit false blocking, a *unit* disagreement is
    error-level but a same-unit *non-overlapping bbox* is only warning-level.
    """
    regions = discover_mesh_regions(case_root)
    if not regions:
        return ()

    diagnostics: list[MeshDiagnostic] = []
    parsed: dict[str, tuple[BoundingBox, ScaleClass]] = {}

    for region in regions:
        label = region.name or "(default)"
        try:
            bbox = read_bounding_box(region.points_path)
        except (MeshParseError, OSError) as exc:
            diagnostics.append(MeshDiagnostic(
                "warning", "mesh_scale_not_checked",
                f"Could not parse points for region {label}: {exc}. "
                f"Run 'checkMeshGeometry' manually to verify scale "
                f"(detect-only by default).",
                region.name,
            ))
            continue
        scale = classify_scale(bbox.max_dim)
        parsed[region.name] = (bbox, scale)
        if scale.unit != "m":
            diagnostics.append(MeshDiagnostic(
                "error", "mesh_not_si",
                f"Region {label} max dimension {bbox.max_dim:g} suggests "
                f"{scale.unit}, not metres. Rescale with "
                f"'checkMeshGeometry -region {region.name or 'region0'} "
                f"-scale {scale.scale_factor:g}'.",
                region.name,
            ))
        elif bbox.max_dim >= _AMBIGUOUS_LOWER:
            diagnostics.append(MeshDiagnostic(
                "warning", "mesh_scale_ambiguous",
                f"Region {label} max dimension {bbox.max_dim:g} m is unusually "
                f"large for a single mesh. Treated as metres (no rescale); "
                f"verify it is an intended large SI domain (e.g. whole torso) "
                f"and not a unit error.",
                region.name,
            ))

    groups = coupled_groups or [frozenset(parsed)]
    for group in groups:
        members = [r for r in group if r in parsed]
        if len(members) < 2:
            continue
        ref = members[0]
        ref_bbox, ref_scale = parsed[ref]
        for other in members[1:]:
            other_bbox, other_scale = parsed[other]
            if other_scale.unit != ref_scale.unit:
                diagnostics.append(MeshDiagnostic(
                    "error", "mesh_scale_mismatch",
                    f"Coupled regions {ref or '(default)'} ({ref_scale.unit}) "
                    f"and {other or '(default)'} ({other_scale.unit}) are at "
                    f"different scales.",
                    other,
                ))
            elif not _bbox_overlaps(ref_bbox, other_bbox):
                diagnostics.append(MeshDiagnostic(
                    "warning", "mesh_scale_mismatch",
                    f"Regions {ref or '(default)'} and {other or '(default)'} "
                    f"have non-overlapping bounding boxes; if they are coupled "
                    f"they may not share a coordinate frame (advisory: coupling "
                    f"is assumed, not derived from the case's own declarations).",
                    other,
                ))

    return tuple(diagnostics)
