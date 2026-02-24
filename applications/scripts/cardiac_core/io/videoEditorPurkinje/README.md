# Purkinje Frame Rendering (ParaView)

Frame rendering is handled by `render_purkinje_frames.py` in this folder. It
renders LV/RV VTU sequences to PNGs and optionally creates GIFs via `ffmpeg`.

## Requirements

- ParaView `pvpython` in your PATH
- `ffmpeg` (only if `--make-gif` is used)

## Quick Start

Run from this folder:

```bash
pvpython render_purkinje_frames.py --make-gif --fps 10 \
  --combined-geometry-vtu ../outputs/<input_stem>_geometry.vtu
```

Defaults:

- LV frames: `frames_biv-line-lv/biv-line-lv_step_*.vtu`
- RV frames: `frames_biv-line-rv/biv-line-rv_step_*.vtu`
- Transparent PNGs, tubes enabled, white background
- Combined output name: `purkinje_generation_geometry`

## FPS Options

- One value applies to LV, RV, and combined:

```bash
--fps 10
```

- Two values use the first for LV (and combined), second for RV:

```bash
--fps 12 8
```

## Outputs

PNG sequences are written under `../outputs/videoEditorPurkinje/png_frames_generation/`:

- `frames_lv_fractal_tree_generation_png/`
- `frames_rv_fractal_tree_generation_png/`
- `frames_purkinje_generation_geometry_png/`

GIFs (when `--make-gif`) are written under `../outputs/videoEditorPurkinje/`:

- `lv_fractal_tree_generation.gif`
- `rv_fractal_tree_generation.gif`
- `purkinje_generation_geometry.gif`

## Geometry Overlay

Geometry is only applied to the combined LV+RV sequence. Pass the geometry VTU with:

```bash
--combined-geometry-vtu ../outputs/<input_stem>_geometry.vtu
```
