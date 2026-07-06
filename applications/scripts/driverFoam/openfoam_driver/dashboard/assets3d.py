from __future__ import annotations

import hashlib
from pathlib import Path

FINAL_3D_PREFIXES = ("PATHOS/", "heartSim3D-1D/")

# Filenames preferred as the primary displayable surface, in priority order.
_SURFACE_PATTERNS = ("*heart*.vtu", "*heart*.vtk", "*.vtu", "*.vtk")


def is_final_3d(case_id: str) -> bool:
    return any(case_id.startswith(p) for p in FINAL_3D_PREFIXES)


def find_primary_surface(case_dir: Path) -> Path | None:
    case_dir = Path(case_dir)
    for pattern in _SURFACE_PATTERNS:
        hits = sorted(case_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def export_glb(surface: Path, out_glb: Path, scalar: str | None = None) -> Path | None:
    """Convert a VTK/VTU/VTP surface to GLB. Returns None if pyvista/trimesh
    are unavailable. Caches on source content hash."""
    surface, out_glb = Path(surface), Path(out_glb)
    hash_file = out_glb.with_suffix(out_glb.suffix + ".hash")
    src_hash = _hash(surface)
    if out_glb.is_file() and hash_file.is_file() and hash_file.read_text() == src_hash:
        return out_glb
    try:
        import pyvista as pv
        import trimesh
    except Exception:
        return None
    out_glb.parent.mkdir(parents=True, exist_ok=True)
    mesh = pv.read(surface)
    surf = mesh.extract_surface() if hasattr(mesh, "extract_surface") else mesh
    surf = surf.triangulate()
    faces = surf.faces.reshape(-1, 4)[:, 1:]
    tm = trimesh.Trimesh(vertices=surf.points, faces=faces, process=False)
    tm.export(out_glb, file_type="glb")
    hash_file.write_text(src_hash)
    return out_glb
