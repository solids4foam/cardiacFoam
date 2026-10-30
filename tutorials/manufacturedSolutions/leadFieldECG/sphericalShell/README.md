# leadFieldECG/sphericalShell

Lead-field ECG against an exact torso potential on a curved heart: a heart
shell `0.5 < |x| < 1` between the bath ball `bloodPool` (`|x| < 0.5`) and the
outer bath `1 < |x| < 1.5`. Current generated in the heart crosses
both spheres into the torso, the heart conductivities are rotated full
tensors whose fibre angle varies with the radius, and the bath conductivity
varies in space. The parent case [leadFieldECG](../README.md) has the
lead-field method and its slab and layered references.

## Problem

The exact fields are the transmission solution of
[`shell_mms.py`](../../generators/shell_mms.py)
(`bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H`):

```text
heart:   div((Gi + Ge) grad phi) = -div(Gi grad Vm) + S_H
torso:   div(sigmaB grad phi)    = S_B
spheres: phi and n.Ge grad phi = n.sigmaB grad phi continuous,
         n.Gi grad(Vm + phi) = 0
outer:   n.sigmaB grad phi = g, the exact current on |x| = 1.5
```

`S_H`, `S_B` and `g` are the solution's own sources. The lead field `z_e` of
electrode `e` solves the transpose of the discrete forward operator, so the
electrode potential is `z_e` applied to the whole right-hand side:

```text
phi_e = sum_cells z_e V (b(Vm) + s) + c_e
s     = S in every cell, minus g |S_f|/V in the owner cell of each outer face
```

with the volume mean of `s` removed (the outer boundary carries no data in
the conductor, so the problem is pure Neumann with a pin). The verifier adds
`s` to the direct forward solve and `sum z_e V s` to the lead-field channel,
for the exact `Vm` at the heart cell centres, and compares both with the
exact potential at the electrodes (rows `leadFieldStatic-*`,
`directStatic-*`). The pinned constant is that of the conductor's
`phiERefPoint`.

The verifier stops when a cell centre lies in another zone of the solution
than the conductor's myocardium/torso split, when Gi, Ge or the torso
conductivity differ from the solution's at the cell centres, or when an
electrode lies in the myocardium.

## Stack

- myocardium: `bidomainSolver` on the cellZone `myocardium`,
  `sealedHeartBoundary true`, `conductivitySource field`, one time step of
  `BuenoOrovio` (the static rows do not use the computed `Vm`);
- conductivities: `0/ConductivityIntracellular`,
  `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity` written by
  `setManufacturedFields` from the solution at the cell centres;
- ECG domain: `leadFieldECG`, torso cellZones `bloodPool` and `bath`,
  `distanceWeightedHarmonic`, pinned in the outer bath; electrodes
  `outerFar`, `outerNear`, `outerSide` (outer bath) and `pool` (bath ball);
- verifier: `manufacturedLeadFieldECGVerifier`, `reference sphericalShell`;
- schemes: `sealedWallLeastSquares` for `Vm` and `phiE`, `leastSquares`
  elsewhere, `Gauss linear corrected`.

## Mesh

- `-mesh hex` (`system/blockMeshDict`): cubed sphere, an inner cube and three
  six-block layers with the faces at `r = 0.5`, `1` and `1.5` projected onto
  the spheres; `N`, a multiple of 4, gives `N` radial heart cells per unit
  length. The cells next to the spheres are nearly orthogonal to them.
- `-mesh tet` (`shell.geo.template`): gmsh Delaunay tetrahedra of size
  `1/N`, interfaces on the spheres; the cells next to the spheres are skewed.

## Run

```bash
./Allrun [parallel] [-mesh hex|tet] [-N n] [-seed n]
```

| option | sets | default |
|---|---|---|
| `parallel` | `decomposePar`, `buildLeadFields` and `cardiacFoam` in parallel | serial |
| `-mesh hex\|tet` | cubed sphere or gmsh tetrahedra | `hex` |
| `-N n` | cells per unit length | 12 |
| `-seed n` | gmsh seed of `-mesh tet` | 1 |

The verifier writes `postProcessing/manufacturedLeadFieldECG_ECG.dat` and
`postProcessing/manufacturedLeadFieldECGSummary_ECG.dat`.

## Studies

`setup/studies/` holds one study per question. Each study is run in copies of
the case, with the commands in its README. `setup/studies/summarize.py` prints
its tables.

| study | question |
|---|---|
| [curvedConvergence](setup/studies/curvedConvergence/README.md) | order of the lead field and of the direct solve against the exact transmission potential at the electrodes and over the field, on the cubed sphere and on tetrahedra |

## Regression

`regression/regressionTest.sh` runs the case defaults serially and compares
the reciprocity gap and the static errors with
`regression/sphericalShell.reference`.
