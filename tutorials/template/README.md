# cardiacFoam Case Template

This directory is a skeleton for creating new cardiacFoam tutorials from scratch.
Copy the whole `template/` tree into your new case directory and fill in the
values described below.

## Directory structure

```
constant/
  physicsProperties     ← always required (1 line: type electroModel)
  electroProperties     ← solver config, ionic model, stimulus, ECG
  [0/]                  ← only needed for pre-existing scalar fields (see below)

system/
  controlDict           ← time control, output frequency
  fvSchemes             ← discretisation schemes (choose A / B / C variant)
  fvSolution            ← linear solver settings (choose A / B / C variant)
  decomposeParDict      ← parallel decomposition (set numberOfSubdomains)
  [blockMeshDict]       ← mesh generator — copy from nearest tutorial and adapt
```

## Required vs optional files

| File | Required? | Notes |
|------|-----------|-------|
| `constant/physicsProperties` | **Always** | `type electroModel;` |
| `constant/electroProperties` | **Always** | Solver, ionic model, stimulus |
| `system/controlDict` | **Always** | Time stepping and output |
| `system/fvSchemes` | **Always** (PDE solvers) | Not needed for singleCellSolver |
| `system/fvSolution` | **Always** (PDE solvers) | Not needed for singleCellSolver |
| `system/decomposeParDict` | Parallel runs | Required if using mpirun |
| `0/` directory | **Only sometimes** | See below |

### When is a `0/` directory needed?

The solver initialises Vm, psi, and all ionic state variables internally —
no `0/Vm` or similar files are required for standard monodomain, bidomain,
eikonal, or singleCell runs.

A `0/` directory **is** needed only when the solver must read a pre-existing
scalar field at startup, for example:

- `0/t` — transmural normalised distance field (0 = endo, 1 = epi), required
  when `ionicHeterogeneity` is configured. See `tutorials/NiedererEtAl2011/tissueNiedererEtAl2011/0/t`.
- Custom fibre orientation or regional fields for advanced setups.

## Mesh

There is no mesh in this template — a mesh must be provided either by:

1. Running `blockMesh` with a `system/blockMeshDict` (for simple slab geometries).
   Copy from `tutorials/NiedererEtAl2011/NiedererEtAl2011verification/system/blockMeshDict`
   and scale `scale`, `vertices`, `blocks`, and `boundary` for your geometry.

2. Importing a pre-existing mesh (VTK, GMSH, etc.) with the appropriate converter.

### Mesh resolution guidelines

| Solver | Recommended cell size | Reason |
|--------|----------------------|--------|
| monodomain explicit | ≤ 0.5 mm | Diffusion stability (CFL) |
| monodomain implicit | ≤ 1.0 mm | Accuracy |
| eikonalSolver | ≤ 0.5 mm | Activation time accuracy |
| bidomainSolver | ≤ 0.5 mm | Coupled Vm + phiE accuracy |

Example: a 20 × 3 × 7 mm slab at 0.5 mm → 40 × 6 × 14 = 3360 cells (fast smoke test).
At 0.1 mm → 200 × 30 × 70 = 420 000 cells (production quality).

## Electrode positions (pseudoECG / eikonalECG)

Electrode positions are **absolute simulation coordinates in metres**.

When using `blockMesh` with `scale 0.001`, a blockMesh vertex at `(20, 3, 7)`
maps to the physical point `(0.02, 0.003, 0.007)` m. Electrodes should be
placed **outside the tissue**, typically 50–200 mm away:

```c++
electrodePositions
{
    lead_I    ( 0.1   0.0   0.0 );   // 100 mm from origin in X
    lead_II   ( 0.0   0.1   0.0 );   // 100 mm from origin in Y
    lead_III  (-0.05  0.05  0.0 );
}
```

For a 20 × 3 × 7 mm slab centred near the origin, placing electrodes at
0.05–0.1 m distance gives a usable pseudo-ECG signal.

## c0 — eikonal conduction velocity

`c0` has units **m/s** (`[0 1 -1 0 0 0 0]`). It is the reference conduction
velocity entering the eikonal wave-speed normalisation.

Physiological range for human ventricular tissue:

| Direction | Typical range |
|-----------|--------------|
| Longitudinal (fibre) | 0.6 – 0.8 m/s |
| Transverse | 0.2 – 0.4 m/s |
| Isotropic approximation | 0.6 – 0.8 m/s |

For a homogeneous slab, `c0 = 0.7 m/s` is a reasonable starting point.
If using an anisotropic conductivity tensor, c0 should reflect the dominant
(fibre-direction) velocity.

## Utilities

All cardiacFoam utilities are catalogued in
`applications/scripts/driverFoam/openfoam_driver/utility_catalog.py`.
Run any of them with `-case <path>` to operate on a specific case directory.

### Mesh preparation

| Utility | When to use |
|---------|-------------|
| `newVtkUnstructuredToFoam <file.vtk>` | Import a 3D VTK unstructured mesh. Also reads scalar/tensor CELL_DATA fields. Use `-no-fields` for mesh-only import. See full workflow below. |
| `transformPoints -scale '(0.001 0.001 0.001)'` *(standard OpenFOAM)* | **Always run after VTK import.** VTK meshes from GMSH/SimNIBS/Meshalyzer are in mm; OpenFOAM needs m. Not in the cardiacFoam catalog. |
| `checkMeshGeometry` | Reports bounding box and can auto-scale units. Run after `transformPoints` to confirm dimensions look correct. |
| `checkMesh` *(standard OpenFOAM)* | Mesh topology/quality: non-orthogonality, skewness, max cell volume ratio. Run after scaling. Not in the cardiacFoam catalog. |
| `1DgraphToFoam <graph.vtk>` | Convert a VTK 1D line network to `constant/purkinjeGraph`. Required before any `conductionNetworkDomains` setup. |

#### Full VTK import workflow

```bash
# 1. Import mesh + fields
newVtkUnstructuredToFoam myHeart.vtk -case ./myCase

# 2. Scale from mm to m (ALWAYS required — VTK has no unit information)
transformPoints -scale '(0.001 0.001 0.001)' -case ./myCase

# 3. Fix Diffusivity dimensions (ALL fields come out as [0 0 0 0 0 0 0])
#    Only Diffusivity needs fixing — fiber, sheet, tags, uvc_transmural are
#    genuinely dimensionless and need no change.
sed -i.bak 's/dimensions.*\[0 0 0 0 0 0 0\]/dimensions      [-1 -3 3 0 0 2 0]/' \
    myCase/0/Diffusivity

# 4. Validate
checkMeshGeometry -case ./myCase   # bounding box sanity
checkMesh -case ./myCase           # topology / quality
```

**On Diffusivity dimensions:** `[-1 -3 3 0 0 2 0]` encodes S/m (conductivity).
If your VTK file exported conductivity in S/mm rather than S/m, multiply the
values by 1000 in addition to fixing the header. Most pipelines use S/m.

### Field setup

| Utility | When to use |
|---------|-------------|
| `setFibreField` | Compute transmural distance `t` and fibre orientation fields (f0, et, en, el) from Laplace solve. Required before `ionicHeterogeneity`. Currently hard-coded for ellipsoidal ventricle geometry (alphaEndo = 60°, alphaEpi = −60°). |
| `setCardiacScarSeverity` | Compute scar depth/severity scalar fields and optionally scale `Diffusivity` in scar cells. Requires a `constant/polyMesh/sets/scarSet` cellSet and `system/cardiacCoreDict/scarSeverity`. |
| `setTorsoOrganConductivityField` | Assign conductivity values per cellZone for torso/bath ECG domains. Produces `0/bodyAndOrgansConductivity` for use with `bathPotentialDomain.bathConductivityField`. |

### Verification (run before full simulation)

| Utility | When to use |
|---------|-------------|
| `listCellModelsVariables` | **Run this first** to see every state, algebraic, constant, and rate variable available for the configured ionic model. Use the output to fill in `outputVariables.ionic.export`. |
| `ionicHeterogeneityProbe` | Validates transmural AP morphology across 0–1 wall depth for BuenoOrovio heterogeneity setups. Outputs APD30/50/70/90, waveform smoothness, and APD envelope reports. Run this before any `ionicHeterogeneity` 3D run. |
| `sweepCurrents` | Plots individual ionic current vs. Vm to validate ionic model setup. Useful when `ionicConstantOverrides` are active — verify the modified current still has physiological shape. Reads `constant/sweepCurrents` config. |
| `runPurkinjeGraph` | Advance the Purkinje graph in isolation without myocardium. Use to validate graph connectivity and initial conditions before coupling with the 3D solver. |

### Post-processing

| Utility | When to use |
|---------|-------------|
| `recomputePseudoECG` | Recompute pseudo-ECG from stored `Vm` fields without rerunning the 3D solver. Use to try new electrode positions or fix a wrong `electrodePositions` block post-hoc. |

### Gap: no electrode-distance utility

There is currently no dedicated utility to check that electrodes are outside the
tissue domain or to report the minimum distance between each electrode and the
nearest heart cell. As a workaround, use `checkMeshGeometry` to read the mesh
bounding box and manually verify that electrode coordinates in
`electrodePositions` are outside that box plus a comfortable margin (≥ 20 mm).

## Allrun pattern

A minimal `Allrun` for a blockMesh-based case:

```bash
#!/bin/bash
cd "$(dirname "$0")"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

runApplication blockMesh
runApplication cardiacFoam
```

For parallel:

```bash
runApplication blockMesh
runApplication decomposePar
runParallel cardiacFoam
runApplication reconstructPar
```

## Checklist before running

- [ ] `constant/physicsProperties` present (`type electroModel;`)
- [ ] `constant/electroProperties` — solver type, ionic model, tissue, chi, cm, conductivity, stimulus
- [ ] `system/controlDict` — endTime covers full AP duration; deltaT matches stability requirement
- [ ] `system/fvSchemes` — correct variant (A = monodomain, B = eikonal, C = bidomain)
- [ ] `system/fvSolution` — correct variant
- [ ] Mesh present (`constant/polyMesh/` or generate via blockMesh)
- [ ] If using ionicHeterogeneity: `0/t` transmural field present
- [ ] If using eikonalECG: `sampling.end` ≥ AP duration (~0.3–0.5 s for ventricle)
- [ ] If parallel: `system/decomposeParDict` present and `numberOfSubdomains` matches mpirun `-np`
