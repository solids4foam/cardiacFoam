# idealizedTorso

The idealized ventricles inside an idealized torso. A heart-only bidomain
source, paced through the Purkinje network, is recorded at the body surface
through lead fields. The case runs the source with either wall treatment of the
sealed heart, so that their effect on the activation, on `phiE` and on the
body-surface ECG can be compared on one geometry.

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
holds them as `(fibre crossFibre)` pairs; `Allrun` runs
`setup/heartConductivity.py`, which writes `0/ConductivityIntracellular` and
`0/ConductivityExtracellular` as `sigmaCross I + (sigmaFibre - sigmaCross) f f`
from the fibre field `0/fiber` (the frame fiber, sheet, sheetNormal is
orthonormal), so the conductivity tensor fields of `../mesh/heartTorso` are not
used. Membrane: `chi` 1400 cm^-1, `cm` 1 µF/cm^2.

## Model

- Source: `bidomainSolver` on the `myocardium` zone, `implicit`, Godunov
  coupling, `backward` in time, `deltaT` 2e-5 s. Bueno-Orovio with the
  transmural regions of `t` (endocardial, M, epicardial) and the apicobasal
  scaling of `tauSi`, as in `../electroHeart`. Purkinje network
  (`monodomain1DSolver`, `purkinjeGraph`), paced at its root, coupled
  one-way to the myocardium.
- Torso: passive, `torsoCoupling leadField`; it does not feed back into the
  source.
- ECG: `torsoECG` at nine body-surface electrodes, E1-E6 (precordial) and RA,
  LA, LL. The precordial leads are V*i* = E*i* - (RA + LA + LL)/3.

## Wall treatment

| Option | `sealedHeartBoundary` | `grad(Vm)`, `grad(phiE)` |
|---|---|---|
| `sealed` (default) | `true` | `sealedWallLeastSquares` |
| `stock` | `false` | `leastSquares` |

`Allrun` writes the choice to `local/wall` and `local/gradSchemes`. The lead
fields depend on it (it enters their signature), so each run builds its own.

## Usage

```bash
./Allrun sealed parallel            # 6 ranks (system/decomposeParDict)
./Allclean
./Allrun stock parallel
```

| Option | Meaning | Default |
|---|---|---|
| `sealed`, `stock` | wall treatment | `sealed` |
| `parallel` | decomposePar, buildLeadFields, cardiacFoam, reconstructPar | serial |
| `-endTime t` | end time (s) | 0.1 (one QRS; the ventricles activate over about 70 ms) |

`buildLeadFields` solves the lead fields through PETSc: cardiacFoam must be
built with solids4foam and PETSc. Outputs: the lead fields and their build
record in `constant/` (or `processor*/constant/`), the source fields every
5 ms and the electrode potentials in `postProcessing/`.

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
