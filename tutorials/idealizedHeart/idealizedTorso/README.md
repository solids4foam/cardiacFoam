# idealizedTorso

The idealized ventricles inside an idealized torso, with the body-surface ECG at
nine electrodes over one beat (QRS and T wave). The bidomain myocardium is
paced through the Purkinje network. The torso is coupled in one of two ways,
chosen at run time:

| Coupling | What is solved | ECG |
|---|---|---|
| `leadField` (default) | heart-only bidomain; the torso does not feed back | lead fields of the heart-torso conductor (`buildLeadFields`), applied to the heart source of `Vm` every step |
| `bath` | one `phiE` over the heart and the torso, solved with `Vm` every step (two-way) | `phiE` of the electrode cells |

Both couplings use the same heart, mesh, ionic regions, Purkinje network, time
step, wall treatment, electrodes and reference point. Their difference is the
effect of the torso on the heart, together with the wall rule below.

## Geometry

`../mesh/heartTorso` (Git LFS) holds the mesh and its `0/` fields. The
heart cells are those of the idealized ventricles of `../mesh`; the heart/torso
interface is conformal.

| Cell zone | Cells | Conductivity (S/m) |
|---|---|---|
| `myocardium` | 123617 | intracellular 0.17 / 0.019, extracellular 0.62 / 0.24 (fibre / cross-fibre) |
| `blood` | 31192 | 0.7 |
| `lungs` | 9680 | 0.0389 |
| `chestWall` | 17027 | 0.05 |
| `torso` | 117966 | 0.22 |

The five zones partition 299482 tetrahedra; the outer patch is `torsoSurface`.
The torso conductivities are the cell field `bodyAndOrgansConductivity`.

The heart conductivities are those of Clerc (1976), as used by Niederer et al.
(2011), transversely isotropic about the fibre. `constant/heartConductivity`
holds them as `(fibre crossFibre)` pairs. `Allrun` runs
`setup/heartConductivity.py`, which writes `0/ConductivityIntracellular` and
`0/ConductivityExtracellular` as `sigmaCross I + (sigmaFibre - sigmaCross) f f`
from the fibre field `0/fiber` (the frame fiber, sheet, sheetNormal is
orthonormal). The conductivity tensor fields of `../mesh/heartTorso` are not
used. Membrane: `chi` 1400 cm^-1, `cm` 1 µF/cm^2.

## Model

- Myocardium: `bidomainSolver` on the `myocardium` zone, `implicit`, Godunov
  coupling, `backward` in time, `deltaT` 2e-5 s.
  - Bueno-Orovio with two transmural regions of `t` (`t` = 0 at the
    endocardium): endocardial on `[0, 0.6]`, epicardial on `[0.6, 1]`, blended
    over a width of 0.1.
  - The apicobasal scaling of `tauSi` (0.97 to 1.04), as in `../electroHeart`.
  - `constant/prePacingProperties` starts each region from its single-cell
    state paced to steady state at 1000 ms.
- Purkinje network: `monodomain1DSolver`, `purkinjeGraph`, paced at its root
  and coupled one-way to the myocardium.
- Torso: passive, no ground and no surface current, so `torsoSurface` is
  insulated. The potential is pinned at the torso `phiERefPoint` with both
  couplings. The heart phiE takes `bathHeartPhiETrace zeroGradient` at the
  heart wall.
- ECG: `torsoECG` at nine body-surface electrodes, E1-E6 (precordial) and RA,
  LA, LL, written every time step. The precordial leads are
  V*i* = E*i* - (RA + LA + LL)/3.

## The heart wall in the two couplings

Given `Vm`, both couplings solve the same extracellular problem:

- `div((Gi + Ge) grad phi) = -div(Gi grad Vm)` in the heart;
- `div(sigma grad phi) = 0` in the torso;
- `phi` and the total current continuous at the heart wall.

They write the wall condition differently, because only the bath solves `Vm`
and `phi` together.

- **`bath`.** The pair (`Vm`, `phiE`) satisfies intracellular insulation,
  `n.Gi grad(Vm + phiE) = 0`. So the current crossing the wall is extracellular
  only, `n.Ge grad(phiE)`. The heart side of the heart-torso faces carries `Ge`.
  The `Vm` equation and the heart block of the `phiE` operator treat the
  intracellular wall faces alike (`insulatedFaceConductivity`; closed on a
  sealed wall).
- **`leadField`.** `Vm` is given and does not satisfy that identity with the
  torso's `phi`. The wall condition is the total current,
  `n.(Gi + Ge) grad(phi) + n.Gi grad(Vm) = n.sigma grad(phi)`. The heart side of
  the faces carries `Gi + Ge`. The wall current `q = n.Gi grad(Vm)` of each wall
  face is moved from the heart cell to the torso cell
  (`passiveConductor::transferWallLayer`). On a sealed heart `q` is small but
  not zero.

For the bath's own `Vm` the two conditions are the same equation, so the lead
fields applied to a bath `Vm` reproduce the bath's `phiE` up to discretisation
error. Using `Ge` with a given `Vm` would drop `n.Gi grad(phi)` at the wall: a
first-order error wherever current crosses it
([`../../manufacturedSolutions/leadFieldECG`](../../manufacturedSolutions/leadFieldECG/README.md),
layered and `wallCurrent` problems).

## Wall treatment

| Option | `sealedHeartBoundary` | `grad(Vm)`, `grad(phiE)` |
|---|---|---|
| `sealed` (default) | `true` | `sealedWallLeastSquares` |
| `stock` | `false` | `leastSquares` |

`Allrun` writes the choice to `local/wall` and `local/gradSchemes`. The lead
fields depend on it (it enters their signature), so each `leadField` run
builds its own.

## Time and output

`endTime` 0.6 s covers the QRS and the T wave.

- With these conductivities half of the myocardium is activated at 0.1 s.
  Activation then grows by 0.62 % of the myocardium per ms, so the last cells
  activate at about 180 ms, and by 200 ms at the latest.
- A single Bueno-Orovio cell paced at 1000 ms repolarises to -80 mV within
  about 290 ms of its stimulus (endocardial, the longer one).
- The last cell has therefore repolarised by about 0.5 s; 0.6 s adds a margin.

The fields are written every 5 ms (121 writes) and the electrode potentials
every time step (`postProcessing/torsoECG.dat`).

## Usage

```bash
./Allrun leadField parallel         # default coupling, sealed, 6 ranks
./Allclean
./Allrun bath parallel
./Allrun bath stock parallel -endTime 0.02
```

| Option | Meaning | Default |
|---|---|---|
| `leadField`, `bath` | torso coupling | `leadField` |
| `sealed`, `stock` | wall treatment | `sealed` |
| `parallel` | `decomposePar`, `buildLeadFields` (`leadField`), `cardiacFoam`, `reconstructPar`, `redistributeRestartState` on the ranks of `system/decomposeParDict` (6, scotch) | serial |
| `-endTime t` | end time (s) | 0.6 |

`Allrun` writes the coupling to `local/coupling` and its myocardium keys to
`local/myocardium`:

- `leadField`: the reference point of the heart-only `phiE`;
- `bath`: `bathPredictorCorrector yes`, so each step solves `Vm` and then the
  two PIMPLE passes (`nOuterCorrectors 2`) of the global `phiE` and `Vm`.

`leadField` needs cardiacFoam built with solids4foam and PETSc
(`buildLeadFields`). `bath` does not.

## Cost

Measured on 6 ranks of one Apple workstation, sealed wall:

| Coupling | Time per step | 0.6 s (30000 steps) | Memory | Output |
|---|---|---|---|---|
| `leadField` | 0.64 s (flat over 0-0.1 s); `buildLeadFields` 41 s | about 5.5 h | below `bath` | about 8 GB |
| `bath` | 1.21 s (flat over 0-0.01 s) | about 10 h | 2.35 GB over 6 ranks | about 11 GB |

In the bath test (`-endTime 0.01`), every global `phiE` solve reached its 1e-11
tolerance in 254-388 iterations. The second pass reduced the `phiE` initial
residual by about 8 times (median 3.3e-3 to 3.4e-4).

## Outputs

- `postProcessing/torsoECG.dat`: the nine electrode potentials every time step.
- Every 5 ms: the heart fields (`Vm`, `phiE`, `activationTime`, ...) and, with
  `bath`, the global `phiE`, `VmGlobal` and `sigmaTotal`. Each write is
  reconstructed, with the ionic restart state.
- With `leadField`: the lead fields and their build record in `constant/` (or
  `processor*/constant/`).

## Folder structure

```text
idealizedTorso/
├── Allrun
├── Allclean
├── README.md
├── constant/
│   ├── electroProperties
│   ├── heartConductivity
│   ├── physicsProperties
│   └── prePacingProperties
├── setup/
│   └── heartConductivity.py
└── system/
    ├── controlDict
    ├── decomposeParDict
    ├── fvSchemes
    └── fvSolution
```
