# NiedererEtAl2011verificationTet — unstructured (tet) slab

Tetrahedral-mesh variant of `../NiedererEtAl2011verification`. It reproduces
the Niederer et al. (2011) monodomain slab benchmark on an **unstructured
tetrahedral** mesh and confirms that it converges to the same total activation
time as the structured-hex case and the published benchmark.

(`NiedererEtAl2012.reference` in the hex case is a regression fixture frozen at
t = 0.04 s, not the published values; this case runs to full activation and is
compared against the hex sweep and the digitized paper data.)

## What differs from the hex case

The geometry (20×3×7 mm slab, fibres along the long x-axis), conductivity,
external stimulus, `tissue epicardialCells`, `fvSchemes`, activation
probes/lines, and `NiedererEtAl2012.reference` are reused from the hex case.
Three things change — the mesh, the ionic model (restored to the benchmark
truth), and the diffusion solve numerics that the tet mesh *requires*:

| | hex case | this case |
|---|---|---|
| mesh generator | `blockMesh` (100×15×35 hexes) | `gmsh` → `gmshToFoam` (isotropic tets) |
| spacing control | cell counts in `blockMeshDict` | characteristic length `lc = dx` in `studies/slab.geo.template` |
| `ionicModel` | `BuenoOrovio` (regression default) | `TNNP` (benchmark reference) |
| `solutionAlgorithm` | `explicit` | `implicit` |
| `nNonOrthogonalCorrectors` | 0 (absent) | `2` |
| Vm linear tolerance | `1e-11` | `1e-15` |

> The hex tutorial ships `ionicModel BuenoOrovio` as a fast regression default;
> the driver overrides it to `TNNP` per-case at runtime (see
> `../NiedererEtAl2011verification/studies/run_manifest.json`, cases
> `explicit_TNNP_epicardialCells_*`). This standalone case bakes in the `TNNP`
> truth so it reproduces the actual Niederer reference, not the regression.

The gmsh box is built directly in **metres** (gmshToFoam applies no scale
factor), folding in the hex case's `blockMesh` `scale 0.001`.

### Why the solve numerics differ (second order on tets)

`fvSchemes` is already tet-safe and is reused unchanged: `gradSchemes
leastSquares` (Gauss-linear gradient is only 1st-order on tets; leastSquares is
2nd-order), `laplacianSchemes Gauss linear corrected`, `snGradSchemes
corrected`. But the `corrected` schemes apply the non-orthogonal correction as
a *deferred* term that must be **iterated** to converge — `nNonOrthogonalCorrectors
≥ 1` — and that iteration only takes effect inside the **implicit** linear
solve. On orthogonal hexes the correction is identically zero, so the hex case
runs `explicit` with 0 correctors and still gets the right answer. On tets that
combination leaves the correction un-converged and the diffusion drops below
second order. This case therefore adopts the same numerics validated to reach
(near) second order on unstructured tets in monodomainPseudoECG's tet overlay (`studies/mesh/tet`, formerly monodomainTetMMS):
`implicit` + `nNonOrthogonalCorrectors 2` + `leastSquares` + tight linear
tolerance. Implicit also removes the explicit-diffusion CFL limit that the fine
tet rungs would otherwise hit at the locked timestep.

## What is being tested

Whether the unstructured (tetrahedral) case reproduces the structured-hex
Niederer result. The timestep is **locked at the coarsest sweep value, 0.01 ms
(`deltaT 1e-5`)**, for every dx, so temporal error stays fixed and only the
mesh refines, over the three rungs:

```
dx = 0.5, 0.2, 0.1 mm   ->   lc = 5e-4, 2e-4, 1e-4 m
```

Each rung runs to full slab activation using the driverFOAM per-dx end times
(`END_TIME_BY_DX` in `openfoam_driver/core/defaults/niederer_2012.py`:
0.5 mm → 0.2 s, 0.2 mm → 0.08 s, 0.1 mm → 0.055 s — coarser meshes have slower
numerical conduction velocity and so need longer to activate the far corner).

**Result: the tetrahedral case converges to the same total activation time as
the structured-hex case and the published Niederer benchmark.** Compare the two
mesh types on *effective* spacing (`dx_eff = (V/nCells)^(1/3)`), not nominal
`lc` — a tet mesh splits each nominal cube into ~5 tets, so its effective
resolution runs ahead of the nominal edge length; a nominal-dx comparison is
misleading.

> Note: the finest rung (dx = 0.1 mm) on this slab is a large tet mesh
> (~millions of cells) and is not needed to show the equivalence — the coarser
> rungs already demonstrate it on matched effective resolution. Start coarse.

## Run modes

Single run (default dx = 0.5 mm, or override `LC`):

```bash
./Allrun                 # serial
./Allrun parallel        # parallel
LC=2e-4 ./Allrun         # dx = 0.2 mm
```

Full dx sweep at the locked timestep:

```bash
./studies/run_tet_sweep.sh                 # dx = 0.5, 0.2, 0.1 mm
DX_VALUES="0.5 0.2" ./studies/run_tet_sweep.sh   # quick two-rung pass
```

Per-dx activation times land in `studies/results/dx_<dx>mm/`.

## Requirements

`gmsh` on `PATH` (same dependency as the monodomainPseudoECG / eikonalECG
tet cases).
