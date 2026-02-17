#!/usr/bin/env pvpython
"""Render LV/RV fractal frames to PNGs and optional GIFs (ParaView pvpython).

Inputs:
- LV/RV VTU frame directories (defaults set in CLI args).
- Optional geometry VTU for combined overlay.

Outputs:
- PNG frames under outputs/videoEditorPurkinje/png_frames_generation/.
- Optional GIFs under outputs/videoEditorPurkinje/.

Main:
- Uses ParaView to render line tubes with configurable camera and colors.
  To add GIFs or change FPS, pass flags when running the script, e.g.:
  pvpython render_purkinje_frames.py --make-gif --fps 10
"""
import argparse
import glob
import subprocess
from pathlib import Path

from paraview.simple import (  # type: ignore
    XMLUnstructuredGridReader,
    GetActiveViewOrCreate,
    GetAnimationScene,
    Hide,
    Render,
    SaveAnimation,
    Show,
    _DisableFirstRenderCameraReset,
)


def parse_vec(value, n=3):
    parts = [float(p.strip()) for p in value.split(",")]
    if len(parts) != n:
        raise ValueError(f"Expected {n} values, got {len(parts)} in '{value}'")
    return parts


def render_sequence(
    file_list,
    color,
    output_dir,
    output_prefix,
    view,
    line_width,
    tubes,
    camera_pos,
    camera_focal,
    camera_up,
    camera_scale,
    transparent,
    geom_reader=None,
    geom_color=None,
    geom_opacity=None,
):
    reader = XMLUnstructuredGridReader(FileName=file_list)
    reader.TimeArray = "None"

    display = Show(reader, view, "UnstructuredGridRepresentation")
    display.Representation = "Surface"
    display.LineWidth = line_width
    display.RenderLinesAsTubes = 1 if tubes else 0
    display.AmbientColor = color
    display.DiffuseColor = color

    if geom_reader is not None:
        geom_display = Show(geom_reader, view, "UnstructuredGridRepresentation")
        geom_display.Representation = "Surface"
        if geom_color is not None:
            geom_display.AmbientColor = geom_color
            geom_display.DiffuseColor = geom_color
        if geom_opacity is not None:
            geom_display.Opacity = geom_opacity

    view.ResetCamera(False, 0.9)
    if camera_pos:
        view.CameraPosition = camera_pos
    if camera_focal:
        view.CameraFocalPoint = camera_focal
    if camera_up:
        view.CameraViewUp = camera_up
    if camera_scale is not None:
        view.CameraParallelScale = camera_scale

    animation_scene = GetAnimationScene()
    animation_scene.UpdateAnimationUsingDataTimeSteps()
    animation_scene.NumberOfFrames = len(file_list)

    Render()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pattern = str(output_dir / f"{output_prefix}.png")
    SaveAnimation(
        output_pattern,
        view,
        ImageResolution=view.ViewSize,
        TransparentBackground=1 if transparent else 0,
    )

    Hide(reader, view)


def render_pair_sequence(
    lv_files,
    rv_files,
    lv_color,
    rv_color,
    output_dir,
    output_prefix,
    view,
    line_width,
    tubes,
    camera_pos,
    camera_focal,
    camera_up,
    camera_scale,
    transparent,
    geom_reader=None,
    geom_color=None,
    geom_opacity=None,
):
    lv_reader = XMLUnstructuredGridReader(FileName=lv_files)
    lv_reader.TimeArray = "None"
    rv_reader = XMLUnstructuredGridReader(FileName=rv_files)
    rv_reader.TimeArray = "None"

    lv_display = Show(lv_reader, view, "UnstructuredGridRepresentation")
    lv_display.Representation = "Surface"
    lv_display.LineWidth = line_width
    lv_display.RenderLinesAsTubes = 1 if tubes else 0
    lv_display.AmbientColor = lv_color
    lv_display.DiffuseColor = lv_color

    rv_display = Show(rv_reader, view, "UnstructuredGridRepresentation")
    rv_display.Representation = "Surface"
    rv_display.LineWidth = line_width
    rv_display.RenderLinesAsTubes = 1 if tubes else 0
    rv_display.AmbientColor = rv_color
    rv_display.DiffuseColor = rv_color

    if geom_reader is not None:
        geom_display = Show(geom_reader, view, "UnstructuredGridRepresentation")
        geom_display.Representation = "Surface"
        if geom_color is not None:
            geom_display.AmbientColor = geom_color
            geom_display.DiffuseColor = geom_color
        if geom_opacity is not None:
            geom_display.Opacity = geom_opacity

    view.ResetCamera(False, 0.9)
    if camera_pos:
        view.CameraPosition = camera_pos
    if camera_focal:
        view.CameraFocalPoint = camera_focal
    if camera_up:
        view.CameraViewUp = camera_up
    if camera_scale is not None:
        view.CameraParallelScale = camera_scale

    animation_scene = GetAnimationScene()
    animation_scene.UpdateAnimationUsingDataTimeSteps()
    animation_scene.NumberOfFrames = max(len(lv_files), len(rv_files))

    Render()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pattern = str(output_dir / f"{output_prefix}.png")
    SaveAnimation(
        output_pattern,
        view,
        ImageResolution=view.ViewSize,
        TransparentBackground=1 if transparent else 0,
    )

    Hide(lv_reader, view)
    Hide(rv_reader, view)


def main():
    parser = argparse.ArgumentParser(
        description="Render LV/RV VTU sequences to PNGs using ParaView."
    )
    parser.add_argument(
        "--lv-dir",
        default="frames_biv-line-lv",
        help="Directory with LV VTU frames.",
    )
    parser.add_argument(
        "--rv-dir",
        default="frames_biv-line-rv",
        help="Directory with RV VTU frames.",
    )
    parser.add_argument(
        "--lv-glob",
        default="biv-line-lv_step_*.vtu",
        help="Glob pattern for LV VTU frames.",
    )
    parser.add_argument(
        "--rv-glob",
        default="biv-line-rv_step_*.vtu",
        help="Glob pattern for RV VTU frames.",
    )
    parser.add_argument(
        "--out-dir",
        default="outputs/videoEditorPurkinje",
        help="Output directory for PNG sequences and GIFs.",
    )
    parser.add_argument(
        "--lv-prefix",
        default="lv_fractal_tree_generation",
        help="Prefix for LV PNG sequence.",
    )
    parser.add_argument(
        "--rv-prefix",
        default="rv_fractal_tree_generation",
        help="Prefix for RV PNG sequence.",
    )
    parser.add_argument(
        "--combined-prefix",
        default="purkinje_generation_geometry",
        help="Prefix for combined LV+RV PNG sequence.",
    )
    parser.add_argument("--width", type=int, default=2196, help="PNG width.")
    parser.add_argument("--height", type=int, default=1240, help="PNG height.")
    parser.add_argument(
        "--lv-color",
        default="1.0,0.3333333333,0.0",
        help="LV RGB color.",
    )
    parser.add_argument(
        "--rv-color",
        default="1.0,1.0,0.0",
        help="RV RGB color.",
    )
    parser.add_argument("--line-width", type=float, default=3.0, help="Line width.")
    parser.add_argument(
        "--no-tubes",
        action="store_true",
        help="Disable rendering lines as tubes.",
    )
    parser.add_argument(
        "--camera-pos",
        default="-30.9516268757,9.049613973,538.2343312713",
        help="Camera position (x,y,z).",
    )
    parser.add_argument(
        "--camera-focal",
        default="31.9712705612,93.9416389465,378.4776153564",
        help="Camera focal point (x,y,z).",
    )
    parser.add_argument(
        "--camera-up",
        default="-0.2693330748,0.8903598877,0.3670408223",
        help="Camera view up (x,y,z).",
    )
    parser.add_argument(
        "--camera-scale",
        type=float,
        default=19.1131521189,
        help="Camera parallel scale.",
    )
    parser.add_argument(
        "--combined-geometry-vtu",
        default=None,
        help="Optional geometry VTU to overlay on the combined LV+RV frames.",
    )
    parser.add_argument(
        "--geometry-color",
        default="0.0,0.0,1.0",
        help="Geometry RGB color.",
    )
    parser.add_argument(
        "--geometry-opacity",
        type=float,
        default=0.2,
        help="Geometry opacity.",
    )
    parser.add_argument(
        "--opaque",
        action="store_true",
        help="Save PNGs with an opaque background.",
    )
    parser.add_argument(
        "--make-gif",
        action="store_true",
        help="Generate GIFs after rendering PNGs (requires ffmpeg).",
    )
    parser.add_argument(
        "--fps",
        type=int,
        nargs="*",
        default=[10],
        help="GIF framerate(s): one value for both, or two values (LV RV).",
    )
    args = parser.parse_args()

    _DisableFirstRenderCameraReset()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames_root = out_dir / "png_frames_generation"
    frames_root.mkdir(parents=True, exist_ok=True)

    lv_files = sorted(glob.glob(str(Path(args.lv_dir) / args.lv_glob)))
    rv_files = sorted(glob.glob(str(Path(args.rv_dir) / args.rv_glob)))
    if not lv_files:
        raise RuntimeError(f"No LV VTU files found in {args.lv_dir}")
    if not rv_files:
        raise RuntimeError(f"No RV VTU files found in {args.rv_dir}")

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [args.width, args.height]
    view.Background = [1.0, 1.0, 1.0]

    camera_pos = parse_vec(args.camera_pos, 3) if args.camera_pos else None
    camera_focal = parse_vec(args.camera_focal, 3) if args.camera_focal else None
    camera_up = parse_vec(args.camera_up, 3) if args.camera_up else None

    geom_reader = None
    geom_color = None
    geom_opacity = None
    if args.combined_geometry_vtu:
        geom_reader = XMLUnstructuredGridReader(FileName=[args.combined_geometry_vtu])
        geom_color = parse_vec(args.geometry_color, 3)
        geom_opacity = args.geometry_opacity

    lv_dir = frames_root / f"frames_{args.lv_prefix}_png"
    rv_dir = frames_root / f"frames_{args.rv_prefix}_png"

    render_sequence(
        lv_files,
        parse_vec(args.lv_color, 3),
        lv_dir,
        args.lv_prefix,
        view,
        args.line_width,
        not args.no_tubes,
        camera_pos,
        camera_focal,
        camera_up,
        args.camera_scale,
        not args.opaque,
    )
    render_sequence(
        rv_files,
        parse_vec(args.rv_color, 3),
        rv_dir,
        args.rv_prefix,
        view,
        args.line_width,
        not args.no_tubes,
        camera_pos,
        camera_focal,
        camera_up,
        args.camera_scale,
        not args.opaque,
    )

    combined_dir = frames_root / f"frames_{args.combined_prefix}_png"
    render_pair_sequence(
        lv_files,
        rv_files,
        parse_vec(args.lv_color, 3),
        parse_vec(args.rv_color, 3),
        combined_dir,
        args.combined_prefix,
        view,
        args.line_width,
        not args.no_tubes,
        camera_pos,
        camera_focal,
        camera_up,
        args.camera_scale,
        not args.opaque,
        geom_reader,
        geom_color,
        geom_opacity,
    )

    if args.make_gif:
        if not args.fps:
            fps_values = [10]
        else:
            fps_values = args.fps
        if len(fps_values) == 1:
            lv_fps = rv_fps = combined_fps = fps_values[0]
        else:
            lv_fps = fps_values[0]
            rv_fps = fps_values[1]
            combined_fps = fps_values[0]
        for prefix, png_dir in (
            (args.lv_prefix, lv_dir),
            (args.rv_prefix, rv_dir),
            (args.combined_prefix, combined_dir),
        ):
            if prefix == args.lv_prefix:
                fps = lv_fps
            elif prefix == args.rv_prefix:
                fps = rv_fps
            else:
                fps = combined_fps
            png_pattern = str(png_dir / f"{prefix}.%04d.png")
            palette = out_dir / f"{prefix}_palette.png"
            gif_path = out_dir / f"{prefix}.gif"
            subprocess.run(
                ["ffmpeg", "-y", "-i", png_pattern, "-vf", "palettegen", str(palette)],
                check=True,
            )
            subprocess.run(
                [
                    "ffmpeg",
                    "-y",
                    "-framerate",
                    str(fps),
                    "-i",
                    png_pattern,
                    "-i",
                    str(palette),
                    "-lavfi",
                    "paletteuse",
                    str(gif_path),
                ],
                check=True,
            )


if __name__ == "__main__":
    main()
