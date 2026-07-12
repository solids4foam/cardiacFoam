# manufacturedSolutions/monodomainTetMMS tutorial

Tetrahedral (unstructured) variant of `monodomainPseudoECG`'s 3D manufactured
monodomain MMS. Reuses the exact solution, diagonal conductivity tensor, and
ionic-model forcing unchanged (identical `constant/` and `system/` dicts);
only the mesh generator changes.

## Purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected` Laplacian
scheme sustains its spatial convergence rate on a genuinely unstructured
tetrahedral mesh of the unit cube:

- a smoothly perturbed hex mesh folds before it reaches high
  non-orthogonality (its ceiling is ~8 deg), whereas tets of the cube reach
  tens of degrees naturally; and
- cardiac anatomical meshes are predominantly tetrahedral, so this is the mesh
  topology the framework must actually cope with.

It closes the "all MMS runs are on orthogonal Cartesian boxes" gap from the
Paper I methods review -- see `[[project_paperI_methods_review]]` memory.

## How the mesh is generated

`setup/box.geo.template` is a gmsh (OpenCASCADE) unit cube with a characteristic
length placeholder `__LC__`. `setup/run_tet_sweep.sh` substitutes `lc = 1/N`
per resolution, meshes with gmsh (legacy msh2 format), and imports via
`gmshToFoam`. All six boundary faces lie on the axis-aligned planes
`x,y,z in {0,1}`, where the manufactured cosine field has zero normal
derivative, so the solver's default zeroGradient (no-flux) boundary stays
compatible with the exact solution -- exactly as on the hex meshes. No `0/`
fields or `createPatch` step is needed: cardiacFoam creates `Vm` with a uniform
zeroGradient boundary (`READ_IF_PRESENT`), so gmsh patch naming is irrelevant.

## Effective mesh spacing and observed order

The manufactured verifier assumes a structured mesh and back-computes an
*effective* spacing `dx = 1/round(cbrt(nCells))` from the total cell count.
For an unstructured tet mesh of the unit cube this is the mean cell size, and
it is the correct convergence abscissa. `setup/summarize_tet.py` therefore
computes the observed order from consecutive `dx` values,
`p = log(e_coarse/e_fine) / log(dx_coarse/dx_fine)`, rather than assuming a
factor-of-two refinement, and reports it next to the `checkMesh` max
non-orthogonality and max skewness so the mesh quality is explicit.

The error norms are the framework's standard cell-count-averaged L2/Linf, used
identically on every MMS case (hex and tet) for consistency; a volume-weighted
norm would be marginally more rigorous on non-uniform tets but is not what the
framework reports.

## Running the sweep

```bash
RESOLUTIONS="10 20 40" ENDTIME=0.02 bash setup/run_tet_sweep.sh
```

- `RESOLUTIONS` -- space-separated nominal cells-per-side (default `10 20 40`).
  `N=80` is available but heavy (~3M tets); add it when you want the fourth
  point.
- `ENDTIME` -- integration window (default `0.02`; the backward `ddt` scheme is
  O(dt^2), so a short window keeps the temporal error subdominant and tet runs
  tractable). Set `ENDTIME=0.2` to match the hex baseline exactly for the paper
  table.

Requires OpenFOAM sourced (the script sources
`/Volumes/OpenFOAM-v2412/etc/bashrc`) and `gmsh` on `PATH`. Results land in
`setup/results/<N>/`, and a combined `setup/results/summary.csv` is written at
the end.
