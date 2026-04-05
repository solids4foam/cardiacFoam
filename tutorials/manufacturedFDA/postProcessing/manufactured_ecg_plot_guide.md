# Manufactured pseudoECG plot guide

This folder contains both monodomain manufactured-solution plots and pseudoECG manufactured-verification plots.

## ECG dimension scope

ECG post-processing is currently generated only for `3D`.

If `1D` or `2D` ECG artifacts are present, they are intentionally ignored with a warning because:
- the numerical pseudoECG in OpenFOAM is accumulated as a 3D cell-volume sum
- the current `1D` and `2D` manufactured ECG references are lower-dimensional integrals

So `1D/2D` remain valid for the manufactured monodomain field verification, but they are not used for ECG verification plots.
During post-processing, archived ECG case files for unsupported dimensions are removed from this folder to keep the ECG artifact set consistent with the 3D-only verification path.

## ECG sweep summary plots

- `manufactured_ecg_electrode_geometry.png`
  Shows the configured 3D electrode locations relative to the unit cube only.

- `manufactured_ecg_electrode_geometry.vtp`
  VTK XML PolyData export of the 3D unit cube and ECG electrode points.
  The file also includes axis lines built from the geometry bounds.
  The point coordinates are normalized to a unit viewing box, and the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  Open this in ParaView if you want to choose your own camera for the geometry view.

- `manufactured_ecg_error_vs_gap_summary.png`
  Two sweep-level 3D-only summary panels versus `N`:
  1. max over electrodes of time-`Linf` numerical/reference error and quadrature difference
  2. mean over electrodes of time-`Linf` numerical/reference error and quadrature difference

  This is the main sweep-level ECG summary plot.
  The numerical/reference error should decrease with refinement, while the quadrature difference should stay much smaller.

- `manufactured_ecg_quadrature_only_summary.png`
  Two quadrature-only sweep panels versus `N`:
  1. max over electrodes of time-`Linf` quadrature difference for each `qCheck`
  2. mean over electrodes of time-`Linf` quadrature difference for each `qCheck`

  This isolates the reference-integral convergence without the numerical pseudoECG error.

- `manufactured_ecg_quadrature_final_time_heatmap.png`
  2D heatmap of the final-time maximum quadrature difference over electrodes:
  - x-axis: `N`
  - y-axis: `qCheck`
  - value: `max_electrodes |reference(qCheck) - reference(qReference)|` at the final time

  This is the direct final-time `N x qCheck` view of the quadrature-only error.

## ECG representative time-series plots

For each dimension/solver pair, the postprocess chooses the highest available `N` case and produces:

- `manufactured_ecg_overlay_<dimension>_<solver>.png`
  Per-electrode time traces of:
  - `numeric`
  - `reference(q=qReference)`

- `manufactured_ecg_quadrature_overlay_<dimension>_<solver>.png`
  Per-electrode time traces of the manufactured reference for all available quadrature orders:
  - `reference(q=6)`
  - `reference(q=12)`
  - `reference(q=24)`
  - `reference(q=48)`
  - `reference(q=qReference)`

  This is the direct “quadratures only” comparison at the electrode points.

- `manufactured_ecg_reference_error_timeseries_<dimension>_<solver>.png`
  Per-electrode semilogy time traces of:
  - `abs error to reference(q=qReference) = |numeric - reference(q=qReference)|`
  - `quadrature difference = |reference(q=qCheck) - reference(q=qReference)|`

  In the raw OpenFOAM columns:
  - `errQ96_<electrode>` means the absolute error to the manufactured reference evaluated with quadrature order 96
  - `deltaQuadratureQ6_Q96_<electrode> = |refQ6_<electrode> - refQ96_<electrode>|`
  - more generally, `errQk_<electrode> = |numeric_<electrode> - refQk_<electrode>|`

- `manufactured_ecg_reference_error_surface_<dimension>_<solver>.png`
  One 3D plot for the representative case with two overlaid colormapped surfaces:
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_primary - refQreference|`

  Surface axes:
  - x-axis: time
  - y-axis: electrode index
  - z-axis: error magnitude
  - z-axis uses a log scale so both surfaces remain visible when their magnitudes differ
  - both surfaces use one shared heatmap/colorbar, and the legend distinguishes them by outline

- `manufactured_ecg_reference_error_surface_<dimension>_<solver>_*.vtp`
  VTK XML PolyData exports of the representative 3D surfaces.
  These files include the surface triangles plus bounding-box and axis line geometry built from the surface bounds.
  The exported point coordinates are normalized to a unit viewing box, while the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  Open these in ParaView (or another VTK viewer) if you want to choose your own camera and render the figure there.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_<electrode>.png`
  Two side-by-side 3D sweep plots for one electrode:
  - left subplot uses the primary check quadrature difference
  - right subplot uses the `q=24` quadrature difference

  In each subplot:
  - x-axis: time
  - y-axis: mesh size `N`
  - z-axis: error magnitude
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_selected - refQreference|`

  The z-axis uses a log scale so the numerical/reference error and the quadrature gap can be read on the same plot without the smaller surface collapsing to zero visually.

  These are the plots to inspect if you want to see how the numerical reference error evolves jointly in time and mesh resolution for a fixed electrode.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_<electrode>_*.vtp`
  VTK XML PolyData exports of the per-electrode sweep surfaces.
  These files include the surface triangles plus bounding-box and axis line geometry built from the sweep bounds.
  The exported point coordinates are normalized to a unit viewing box, while the original coordinates are stored as point-data arrays `raw_x`, `raw_y`, and `raw_z`.
  These use the raw surface values, not the Matplotlib camera or the PNG view.

- `manufactured_ecg_sweep_surface_<dimension>_<solver>_all_electrodes.png`
  Multi-panel overview with all electrodes in one figure.
  Each row is one electrode and contains two 3D plots:
  - left column: primary check quadrature difference
  - right column: `q=24` quadrature difference

  In each subplot:
  - viridis surface with black mesh outline: `|numeric - refQreference|`
  - viridis surface with white mesh outline: `|refQcheck_selected - refQreference|`

## Configured quadrature orders

- `qChecks = 6 12 24 48`
- `qReference = 96`

All ECG plots in this folder are derived only from OpenFOAM outputs built from those configured quadrature orders.
