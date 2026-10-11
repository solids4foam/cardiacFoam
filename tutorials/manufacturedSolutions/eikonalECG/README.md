# eikonalECG

Eikonal-diffusion manufactured solution on the unit cube, on Cartesian
hexahedra (`-mesh hex`, default) or gmsh Delaunay tetrahedra (`-mesh tet`),
with the eikonal-driven ECG checked against its reference integral. The exact
activation time is imposed on every face and on the cells whose centre has
the smallest x, y or z (on hexahedra the first cell layer), so the case checks
the interior discretisation (the metric term, the corrected Laplacian, the
cell-gradient reconstruction, the outer fixed point) and the template ECG; it
does not check a wall treatment, the zero-gradient outer boundary of ordinary
runs or the placement of internal seeds. Sealed curved walls are checked by
[sphericalShell](sphericalShell/README.md).

Origin: this work, for both the activation-time and the eikonal-driven ECG
problems.

## Problem

    -div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1 + S    in [0,1]^d
    psi = psi_ex    on every face and on the cells of smallest centre x, y or z

`psi_ex = exp(k . x)` with `k = (1,0,0)`, `(1,1,0)`, `(1,1,1)` in 1D, 2D,
3D, `M = G/(chi Cm)`, `chi = 3`, `Cm = 2`, `c0 = 1`, and
`S = psi_ex (c0 sqrt(K) - K) - 1`, `K = k . M k`. `G` is

- rotated (`-conductivity rotated`, default): the full tensor `(0.106875
  -0.0084931 -0.022561 0.116682 -0.0130256 0.0398782)`
- axis (`-conductivity axis`): `diag(0.111453302, 0.121585420, 0.030396355)`

The ECG is the lead-field integral of the template `Vm = U(t - psi)`,
`U(s) = sin(2 pi s)`, at external electrodes:

    phi(x_e, t) = int G grad Vm . (x - x_e)/|x - x_e|^3 dx,   grad Vm = -U'(t - psi) grad psi

sampled every 0.01 for `t` in `[0, 2]`. The solver evaluates it with the
reconstructed `grad psi` of the computed field; the reference is a
Gauss-Legendre order-48 quadrature of the exact integrand.

### Stack

- myocardium solver: `eikonalSolver`, the metric term `c0 |grad psi|_M` as
  the bounded convection term `u . grad psi` (`div(phiU,activationTime) Gauss
  linearUpwind limitedGrad`) or evaluated explicitly (`-formulation`)
- outer fixed point: `nOuterCorrectors 8000`, relaxation 0.9, residual
  tolerance `1e-8` (`-tolerance`)
- field verifier: `manufacturedEikonalVerifier`
- ECG: `eikonalECG` with `manufacturedEikonalECGVerifier`
- gradient: `leastSquares` or `Gauss linear` (`-gradient`)

## Mesh

- `-mesh hex` (`system/blockMeshDict.<d>D`): `N` cells per unit length along
  each of the `d` directions.
- `-mesh tet` (`box.geo.template`, 3D): Delaunay tetrahedra of size `1/N` from
  `gmsh -setnumber lc 1/N -setnumber seed <seed>`, imported with `gmshToFoam`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/` and the mesh. `Allrun` stops with
`run ./Allclean first` if a mesh log exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (6 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-mesh hex\|tet` | blockMesh box or gmsh tetrahedra (3D) | hex |
| `-dim 1\|2\|3` | dimension of the hex box; the ECG runs in 3D only | 3 |
| `-N n` | cells per unit length; tet: element size `1/n` | 20 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-conductivity rotated\|axis` | conductivity tensor | rotated |
| `-formulation advectionDiffusion\|explicit` | metric term as a bounded convection term or explicit | advectionDiffusion |
| `-gradient leastSquares\|GaussLinear` | cell-gradient reconstruction | leastSquares |
| `-tolerance x` | outer fixed-point residual tolerance | `1e-8` |
| `-electrodes five\|all` | ECG at E1-E5, or at E1-E5 and the 156 random electrodes R1-R156 | five |
| `-errorField` | write `activationTimeError` and the cell centres to the final time | off |

E1-E5 lie outside the cube; R1-R156 are the electrodes of
[monodomainPseudoECG](../monodomainPseudoECG/README.md). The ECG runs in 3D
only.

`postProcessing/manufacturedEikonalActivationTime.dat` holds the
activation-time errors and the bulk, boundary and gradient rows of
[Bulk, boundary and gradient rows](../README.md#bulk-boundary-and-gradient-rows);
in 3D, `postProcessing/eikonalECG.dat` holds the electrode traces,
`postProcessing/manufacturedEikonalECG_ECG.dat` the traces with their
references, and `postProcessing/manufacturedEikonalECGSummary_ECG.dat` the L1,
L2 and Linf error over time of every electrode.

Regression (the defaults: 3D hex, rotated conductivity, `N = 20`, parallel on
six ranks): the activation-time errors with their bulk, boundary and gradient
rows, the ECG errors of E1-E5 and their final trace values.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables.

| study | question |
|---|---|
| [cartesianConvergence](setup/studies/cartesianConvergence/README.md) | spatial order of the activation time on hexahedra, 1D/2D/3D, and of the 3D ECG |
| [tetConvergence](setup/studies/tetConvergence/README.md) | spatial order on tetrahedra: `leastSquares` against `Gauss linear`, axis against rotated conductivity, the two formulations; isolated gradient reconstruction |
| [nonlinearControl](setup/studies/nonlinearControl/README.md) | share of the outer stopping tolerance in the tetrahedral errors |
| [errorLocalisation](setup/studies/errorLocalisation/README.md) | bulk and boundary-adjacent parts of the tetrahedral activation error, both gradients |

## Variants

- [sphericalShell](sphericalShell/README.md): sealed curved walls crossed by a
  rotated conductivity, with a pinned stimulus ball.
