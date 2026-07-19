# Standalone bidomain corrector study

This study separates two operations that were previously both expressed through
PIMPLE iteration counts on tetrahedral meshes:

- an outer bidomain sweep, which repeats the coupled `phiE -> Vm` block; and
- an equation-level non-orthogonal reassembly, which updates the deferred
  corrected-Laplacian contribution within each equation.

The four variants cross `nOuterCorrectors = 1,2` with
`nNonOrthogonalCorrectors = 0,1`. They use the same unit-cube Delaunay mesh
family, least-squares gradient, corrected Laplacian, and `dt ~ h^2` controls as
the monodomain tetrahedral verification. The short fixed-step window is a
same-mesh iteration-sensitivity test, not a replacement for the four-level
Cartesian bidomain convergence study.

Run from the repository root:

```bash
bash tutorials/manufacturedSolutions/bidomain/setup/correctorStudy/run_corrector_study.sh
```

Override `RESOLUTIONS`, `VARIANTS`, `RESULTS_DIR`, or `KEEP_WORK` through the
environment for focused reruns.

