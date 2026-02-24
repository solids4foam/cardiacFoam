#!/usr/bin/env python3
"""Tag endocardial/epicardial surfaces from UVC fields."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "cardiac_preproc" / "src"
DEFAULT_CONFIG = Path(__file__).resolve().with_name("tag_endo_epi_surface_config.py")

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cardiac_preproc.tagging.tag_endo_epi_surface import main  # noqa: E402


def _load_config(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load config from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _extract_input_arg(args: list[str]) -> str | None:
    for i, token in enumerate(args):
        if token == "--input" and i + 1 < len(args):
            return args[i + 1]
        if token.startswith("--input="):
            return token.split("=", 1)[1]
    return None


def _has_flag(args: list[str], name: str) -> bool:
    return any(token == name for token in args)


def _has_seed_arg(args: list[str]) -> bool:
    return any(token == "--seed-point" or token.startswith("--seed-point=") for token in args)


def _pick_seed_point(input_vtk: str) -> str:
    import numpy as np
    import pyvista as pv

    mesh = pv.read(input_vtk).extract_surface().triangulate()
    if mesh.n_cells == 0:
        raise RuntimeError("Input has no surface cells to pick from.")

    points = np.asarray(mesh.points)
    if points.size == 0:
        raise RuntimeError("Input surface has no points to pick from.")
    picked: dict[str, np.ndarray | None] = {"seed": None}

    def on_pick(point):
        if point is None:
            return
        p = np.asarray(point, dtype=float).reshape(-1)
        if p.size != 3:
            return
        idx = int(np.argmin(np.linalg.norm(points - p, axis=1)))
        snapped = np.asarray(points[idx], dtype=float)
        picked["seed"] = snapped
        print(
            "Clicked point: "
            f"({p[0]:.6f}, {p[1]:.6f}, {p[2]:.6f}) | "
            "Nearest mesh point: "
            f"({snapped[0]:.6f}, {snapped[1]:.6f}, {snapped[2]:.6f})"
        )

    pl = pv.Plotter()
    pl.add_text("Click seed point on surface, then press q", font_size=11)
    pl.add_mesh(mesh, color="lightgray", opacity=0.35, show_edges=False)
    pl.enable_surface_point_picking(
        callback=on_pick,
        left_clicking=True,
        show_point=True,
        show_message=True,
    )
    pl.show()

    seed = picked.get("seed")
    if seed is None:
        raise RuntimeError("No seed point selected.")
    return ",".join(f"{float(x):.6f}" for x in seed)


if __name__ == "__main__":
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", default=str(DEFAULT_CONFIG), help="Path to tagging config file.")
    pre_parser.add_argument(
        "--pick-seed",
        action="store_true",
        help="Open an interactive window to pick --seed-point on the surface.",
    )
    pre_args, remaining = pre_parser.parse_known_args()
    cfg = _load_config(Path(pre_args.config))

    if pre_args.pick_seed and not _has_seed_arg(remaining):
        input_vtk = _extract_input_arg(remaining)
        if not input_vtk:
            raise SystemExit("`--pick-seed` requires `--input <mesh.vtk>`.")
        seed_str = _pick_seed_point(input_vtk)
        print(f"Picked seed point: {seed_str}")
        remaining.append(f"--seed-point={seed_str}")

    config_seed = getattr(cfg, "SHARED_BOUNDARY_SEED", None)
    if config_seed and not _has_seed_arg(remaining):
        remaining.append(f"--seed-point={config_seed}")

    if not getattr(cfg, "COMPLEMENTARY_SURFACE_OUTPUT", True) and not _has_flag(remaining, "--no-surface-output"):
        remaining.append("--no-surface-output")
    if not getattr(cfg, "EXTRACT_VOLUME_OUTPUT", False) and not _has_flag(remaining, "--no-volume-output"):
        remaining.append("--no-volume-output")

    sys.argv = [sys.argv[0], *remaining]
    main()
