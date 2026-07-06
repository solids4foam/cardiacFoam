from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from .assets3d import find_primary_surface

_CANDIDATES = [
    "/Applications/Blender.app/Contents/MacOS/Blender",
    "/Applications/Blender 2.app/Contents/MacOS/Blender",
]
_BLENDER_SCRIPT = Path(__file__).parent / "_blender_render.py"


def find_blender(explicit: str | None) -> str:
    import os
    import shutil
    for cand in [explicit, os.environ.get("BLENDER"), *(_CANDIDATES),
                 shutil.which("blender")]:
        if cand and Path(cand).exists():
            return cand
    raise FileNotFoundError(
        "Blender not found. Pass --blender <path> or set $BLENDER."
    )


def _to_stl(surface: Path, out_stl: Path) -> Path:
    """Convert VTK/VTU to STL via pyvista so Blender can import it."""
    import pyvista as pv
    out_stl.parent.mkdir(parents=True, exist_ok=True)
    mesh = pv.read(surface)
    surf = mesh.extract_surface() if hasattr(mesh, "extract_surface") else mesh
    surf.triangulate().save(out_stl)
    return out_stl


def render(case_dir: Path, out_dir: Path, blender: str | None,
           scalar: str | None = None) -> Path:
    case_dir, out_dir = Path(case_dir), Path(out_dir)
    surface = find_primary_surface(case_dir)
    if surface is None:
        raise FileNotFoundError(f"no VTK/VTU surface under {case_dir}")
    blender_bin = find_blender(blender)
    out_dir.mkdir(parents=True, exist_ok=True)
    stl = _to_stl(surface, out_dir / "surface.stl") if surface.suffix != ".stl" else surface
    out_blend = out_dir / "scene.blend"
    out_png = out_dir / "render.png"
    subprocess.run(
        [blender_bin, "--background", "--python", str(_BLENDER_SCRIPT),
         "--", str(stl), str(out_blend), str(out_png)],
        check=True,
    )
    return out_png


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Render a cardiacFoam 3D case in Blender.")
    p.add_argument("--case", required=True, type=Path)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--blender", default=None)
    p.add_argument("--scalar", default=None)
    args = p.parse_args(argv)
    out = args.out or (Path("dashboard/renders") / args.case.name)
    png = render(args.case, out, args.blender, args.scalar)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
