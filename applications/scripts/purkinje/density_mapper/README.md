# Purkinje Density Probability Engine (Skeleton)

Standalone, extensible scaffold to assign a probability to anatomical ventricular locations while minimizing code duplication.

## Design goals

- Reuse existing tagging/geometry functions through a single feature-extraction interface.
- Keep probability logic modular (`tag -> model` registry).
- Keep implementation standalone (pure Python standard library).

## Package layout

- `src/purkinje_density/types.py`: shared dataclasses and type aliases.
- `src/purkinje_density/features.py`: feature extraction utilities.
- `src/purkinje_density/models.py`: model implementations.
- `src/purkinje_density/registry.py`: tag-to-model registry.
- `src/purkinje_density/engine.py`: high-level API.
- `src/purkinje_density/defaults.py`: build engine from config dict.
- `src/purkinje_density/tree_metrics.py`: tree presence + density metrics.
- `src/purkinje_density/interactive_bullseye.py`: interactive segment-% input UI.
- `src/purkinje_density/vtk_bullseye_workbench.py`: VTK + 3D + interactive bullseye mapping.

## Quick start

```python
from purkinje_density import LocationSample
from purkinje_density.defaults import build_engine_from_dict

config = {
    "models": {
        "lv_apical_septal": {
            "type": "linear_logistic",
            "weights": {"uvc_longitudinal": -3.0, "is_lv": 1.5, "is_septal": 1.2},
            "bias": 1.0,
        },
        "rv_mid_freewall": {
            "type": "range_rule",
            "bias": -0.2,
            "rules": [
                {"feature": "uvc_longitudinal", "min": 0.3, "max": 0.7, "weight": 1.0},
                {"feature": "is_rv", "min": 1.0, "max": 1.0, "weight": 1.3},
                {"feature": "is_freewall", "min": 1.0, "max": 1.0, "weight": 1.1},
            ],
        },
    }
}

engine = build_engine_from_dict(config)

sample = LocationSample(
    chamber="LV",
    wall="septal",
    scalar_fields={"uvc_longitudinal": 0.15},
)

p = engine.probability("lv_apical_septal", sample)
print(p)
```

## Minimal Start: Presence Probabilities Only

`compute_tree_presence_metrics(...)` computes:

1. Probability of any tree part in endocardium/myocardium.
2. Probability of tree terminal ends in endocardium/myocardium.

Example: `examples/tree_presence_usage.py`

## Interactive Bullseye Input

`collect_segment_percentages_interactive(...)` opens a clickable bullseye:

1. Click a segment.
2. Enter its `%` in the figure-side input panel.
3. Segment updates live in the bullseye.
4. Close figure to finish (missing segments auto-fill with `0.0` by default).
5. Visualization tweak: apical ring starts at `45°` (from `0°`) to avoid left-line alignment.
6. Orientation labels are shown: left `SEPTAL`, right `LATERAL`, top `POSTERIOR`, bottom `ANTERIOR`.
7. Values are per-segment (no global sum constraint) and default scale is `0..10`.

Example: `examples/interactive_bullseye_input.py`

## VTK + 3D + Bullseye Workbench

Use one command to:

1. Open a 3D view of the VTK anatomical segmentation.
2. Open the interactive bullseye for `%` input per segment.
3. Map `%` values back to the VTK as a new scalar field.
4. Optionally save output VTK.

Command:

```bash
PYTHONPATH=src python3 -m purkinje_density.vtk_bullseye_workbench \
  --vtk /path/to/mesh.vtk \
  --config-source anatomical \
  --reference aha17 \
  --segment-field anatomical_tag \
  --segment-association cell \
  --input-mode widget \
  --value-max 10 \
  --mapped-field bullseye_percentage \
  --output /path/to/mesh_bullseye_mapped.vtk
```

Notes:

- Workbench visualization requires `pyvista` and `matplotlib` (and typically `numpy` via PyVista stack).
- `--config-source anatomical` reads `../anatomicalSegmentationModels/configs/lv_division_config.py` by default.
- If that file is unavailable, it falls back to local `aha17`.
- Segment field is auto-detected if omitted (`anatomical_tag`, `tagged_segment_id`, ...).
- If your VTK is not anatomically tagged yet, run `anatomicalSegmentationModels/anatomical_segmentModels.py` first, then run this workbench.

## Unified GUI (Single Window, 2 Panels)

Use a single window containing:

1. Interactive bullseye editor with three tabs:
   - Purkinje Density
   - PMJ /mm2
   - Thickness Bundles
2. Mapped 3D values preview (active tab field) with segment-id labels at barycenters.

Behavior:

- `Set Value` edits selected segment on bullseye.
- `Apply to 3D` refreshes the 3D preview (only action that updates 3D).
- `Save + Exit` writes all three user-defined fields plus the automatic area field to output VTK and closes the GUI.
- Output also includes an automatic field: `endocardial_surface_area_by_tag`.
- Default ranges: Purkinje `0..1`, PMJ `/mm2` `0..1`, Thickness Bundles `0..200` (independent scales).

Command:

```bash
PYTHONPATH=src python3 -m purkinje_density.unified_workbench_gui \
  --vtk /path/to/mesh.vtk \
  --config-source anatomical \
  --reference apical4_mid5_basal5 \
  --segment-association auto \
  --purkinje-max-value 1 \
  --pmj-max-value 1 \
  --thickness-bundles-max-value 200 \
  --purkinje-field purkinje_density \
  --pmj-field pmj_per_mm2 \
  --thickness-bundles-field thickness_bundles \
  --endocardial-area-field endocardial_surface_area_by_tag \
  --output /path/to/mesh_bullseye_mapped.vtk
```

When `--output` is set, the GUI also writes a sidecar JSON:

- `/path/to/mesh_bullseye_mapped.vtk.solution.json`

You can re-apply that JSON later to any tagged mesh with:

```bash
PYTHONPATH=src python3 -m purkinje_density.apply_solution_json \
  --vtk /path/to/tagged_mesh.vtk \
  --solution-json /path/to/mesh_bullseye_mapped.vtk.solution.json \
  --segment-association auto \
  --purkinje-field purkinje_density \
  --pmj-field pmj_per_mm2 \
  --thickness-bundles-field thickness_bundles \
  --endocardial-area-field endocardial_surface_area_by_tag \
  --output /path/to/tagged_mesh_with_densities.vtk
```

Optional endocardial selector (if your mesh has a dedicated field):

```bash
  --endocardial-selector-field tissue \
  --endocardial-selector-value endocardium
```

If selector is omitted, code auto-detects `endo_surface` (or `endo_surface_lv`/`endo_surface_rv`) and only falls back to all cells when those tags are absent.

Dependencies:

- `numpy`, `pyvista`, `pyvistaqt`, `matplotlib`, and either `PySide6` or `PyQt5`.

## Tree Density Metrics (Optional Later Layer)

`compute_tree_density_metrics(...)` computes:

1. Probability of any tree part in endocardium/myocardium.
2. Probability of tree terminal ends in endocardium/myocardium.
3. Terminal-end density by tissue in `ends/mm^2`.
4. Terminal-end density by tagged section in `ends/mm^2`.

Example: `examples/tree_metrics_usage.py`

## Run tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py'
```
