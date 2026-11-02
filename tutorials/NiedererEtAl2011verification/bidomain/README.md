# NiedererEtAl2011verification/bidomain

The Niederer et al. (2011) slab benchmark solved with the bidomain model: the
geometry, stimulus, ionic model, probe points and probe lines of the
[monodomain case](../README.md), with `bidomainSolver` in place of
`monodomainSolver`.

- Electro model: `bidomainSolver`, `solutionAlgorithm implicit`
- Ionic model: `TNNPcompactBatched` with `batchedIntegrator rushLarsen`,
  `tissue epicardialCells`
- Outputs: `activationTime` at the nine `system/Niedererpoints` probes and along
  the `system/Niedererlines` diagonal

## Problem

```text
chi Cm dVm/dt = div(G_i grad Vm) + div(G_i grad phiE) - chi Cm Iion + I_s
0             = div(G_i grad Vm) + div((G_i + G_e) grad phiE)
n . G_i grad(Vm + phiE) = n . G_e grad phiE = 0             on every wall
```

on the 20 x 3 x 7 mm slab, with fibres along x, `chi = 140000 1/m`,
`Cm = 0.01 F/m^2` and the stimulus of the monodomain case (50000 A/m^3 for
2 ms in the 1.5 x 1.5 x 1.5 mm corner block `x < 1.5 mm`, `y < 1.5 mm`,
`z > 5.5 mm`). `G_i` and `G_e` are diagonal, `diag(l, t, t)`.

### Conductivity sets (`-conductivity`)

`published` is the bidomain conductivity of Niederer et al. (2011):

| tensor | longitudinal | transverse |
|---|---|---|
| `G_i` | 0.17 S/m | 0.019 S/m |
| `G_e` | 0.62 S/m | 0.24 S/m |

The monodomain case uses their harmonic combination
`G_m = G_i (G_i + G_e)^-1 G_e`, `0.17*0.62/0.79 = 0.1334177215` and
`0.019*0.24/0.259 = 0.01760617761` S/m (10 significant digits).

`equal` (the default) has one anisotropy ratio, `G_e = lambda G_i` with
`lambda = 0.62/0.17`, and the monodomain tensor `G_m` of the monodomain case
as its harmonic combination:

```text
G_i* = G_m (0.17 + 0.62)/0.62
G_e* = G_m (0.17 + 0.62)/0.17
G_i* G_e* / (G_i* + G_e*) = G_m
```

| tensor | longitudinal | transverse |
|---|---|---|
| `G_i*` | 0.16999999997580645 S/m | 0.022433677922419355 S/m |
| `G_e*` | 0.6199999999117647 S/m | 0.08181694301117647 S/m |

These are the double-precision values of the formulas with `G_m = (0.1334177215,
0.01760617761)`; their harmonic combination rounds to `G_m` in double
precision, and `G_e*/G_i*` is the same double in both directions:

```bash
python3 -c "
m = (0.1334177215, 0.01760617761)
i = [g*(0.17 + 0.62)/0.62 for g in m]
e = [g*(0.17 + 0.62)/0.17 for g in m]
print(i, e, [a*b/(a + b) for a, b in zip(i, e)], [b/a for a, b in zip(i, e)])"
```

### What the equal set reproduces

With `G_e = lambda G_i` the elliptic equation gives
`phiE = -Vm/(1 + lambda) + c`, and the parabolic equation becomes the
monodomain equation with `G_i lambda/(1 + lambda) = G_m`. The finite-volume
equations keep this identity: the `phiE` matrix is `(1 + lambda)` times the
`G_i` Laplacian, and the reference cell fixes only `c`. The equal set at a
given `-dx`, `-dt` and `-endTime` reproduces the activation times of the
monodomain case on the same mesh and time step to the tolerances of the
linear solves and of the outer `phiE`-`Vm` iterations. The
[benchmark ladder](setup/studies/benchmarkLadder/README.md) gives the
measured differences.

## Folder structure

```text
tutorials/NiedererEtAl2011verification/bidomain/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   ├── fvSolution
│   ├── Niedererlines
│   └── Niedererpoints
├── setup/
│   └── studies/
│       └── benchmarkLadder/
├── regression/
│   ├── regressionTest.sh
│   ├── equal.reference
│   └── published.reference
├── Allrun
├── Allclean
└── README.md
```

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include; with
no options the case runs the dictionary defaults. `Allclean` removes `local/`
and the mesh. `Allrun` stops with `run ./Allclean first` if a mesh log
exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (6 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar`, `redistributeRestartState` | serial |
| `-dx x` | isotropic cell size in metres; `0.02/x`, `0.003/x` and `0.007/x` must be integers (`5e-4`: 40 x 6 x 14, `1e-4`: 200 x 30 x 70 cells) | `2e-4` (100 x 15 x 35 cells) |
| `-dt x` | time step in seconds | `1e-5` |
| `-endTime x` | end time in seconds, also the write interval | `0.015`, written every `5e-3` |
| `-conductivity equal\|published` | conductivity set | `equal` |

The defaults are the mesh and time settings of the monodomain case, so
`./Allrun` and `../Allrun` solve the same discrete problem with the two
models.

After the run, `postProcess -func Niedererpoints -latestTime` and
`-func Niedererlines -latestTime` write
`postProcessing/Niedererpoints/0/activationTime` and
`postProcessing/Niedererlines/0/activationTime`; `-1` is a point not yet
activated.

## Regression

`regression/regressionTest.sh` runs both conductivity sets at `-dx 5e-4`
(40 x 6 x 14 cells), `-dt 1e-5` and `-endTime 0.015`, `equal` in parallel and
`published` in serial, and compares the activation times of probe 0, of
probes 4 and 8 (not yet activated) and of the line points 5, 8 and 10 against
`regression/equal.reference` and `regression/published.reference`. On this
mesh the wavefront has reached line point 14 of 100 at 15 ms; the two sets differ by up to 7e-6 s at line point 10,
and the references hold to 1e-6 s. The regression takes about two minutes.

```bash
./regression/regressionTest.sh
```

## Studies

| study | question |
|---|---|
| [benchmarkLadder](setup/studies/benchmarkLadder/README.md) | the Niederer resolutions Δx = 0.5, 0.2, 0.1 mm at Δt = 0.005 ms for both conductivity sets, and the equal set against the monodomain case |

## Notes on dictionary entries

`constant/electroProperties` holds the equal set; `-conductivity published`
writes `local/conductivity`, included at the end of `bidomainSolverCoeffs`.
`phiERefPoint (0.01995 0.00295 0.00005)` lies in the corner cell farthest from
the stimulus at every `-dx` of the ladder, away from cell faces, and
`phiEReferenceValue 0` fixes `phiE` there. `sealedHeartBoundary false` and the
`leastSquares` gradient are those of the monodomain case: the slab walls are
normal to the axes of the diagonal tensors, so the wall fluxes are zero
without the sealed-wall treatment.

`system/fvSolution` solves `phiE` with PCG/DIC to an absolute tolerance of
`1e-15`, as the bidomain manufactured solutions, with a cap of 5000 PCG
iterations per solve. The `PIMPLE` loop re-solves `phiE` and `Vm` on every
outer corrector until the initial residuals of both are below `1e-6`
(`residualControl`), with at most 1000 outer correctors and two
non-orthogonal correctors, as the monodomain case.

`system/blockMeshDict` reads the cell counts `_nx`, `_ny`, `_nz` from
`local/mesh`; `system/controlDict` reads `deltaT`, `endTime` and
`writeInterval` from `local/time`.
