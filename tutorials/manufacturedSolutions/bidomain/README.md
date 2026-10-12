# bidomain

Bidomain manufactured solution on the unit cube, on Cartesian hexahedra
(`-mesh hex`, default, 1D/2D/3D) or gmsh Delaunay tetrahedra (`-mesh tet`).
The field has a zero normal gradient on every face and the conductivities are
diagonal, so this case checks the interior operators of the one-mesh bidomain
stack (the elliptic `phiE` solve, the parabolic `Vm` solve, their outer
coupling sweeps, cell-gradient reconstruction, non-orthogonal correction,
ionic coupling, time integration), not a wall treatment; sealed curved walls
are checked by [sphericalShell](sphericalShell/README.md).

Origin: the bidomain problem of the FDA regulatory science tool (Pathmanathan
& Gray 2014).

## Problem

```text
chi Cm dVm/dt = div(G_i grad Vm) + div(G_i grad phiE) - chi Cm Iion + I_s
0             = div(G_i grad Vm) + div((G_i + G_e) grad phiE)
n . grad Vm = n . grad phiE = 0                             on every face
```

in `[0,1]^d`, with `phiI = Vm + phiE` and

- `Vm = sqrt(1 + t) F`, `F = cos(pi x)`, `cos(pi x) cos(2 pi y)`,
  `cos(pi x) cos(2 pi y) cos(3 pi z)` in 1D, 2D, 3D
- `phiE = -k Vm + c`, `phiI = (1 - k) Vm + c`, `k = 1/sqrt(2)`; `c` is a
  constant gauge
- `G_i = diag(0.111453302, 0.121585420, 0.030396355)`,
  `G_e = ((1 - k)/k) G_i = diag(0.046165469, 0.050362330, 0.012590582)`, so the
  elliptic equation holds with no `phiE` source

with `chi = 3`, `Cm = 2` and the FDA ionic closure (`u1`, `u2`, `u3`). The
solver fixes the gauge of `phiE` in the cell containing `phiERefPoint (0.125
0.05 0.05)`. The verifier reports `phiE` and `phiI` against the exact field
with `phiE = phiEReferenceValue` at that point, and `phiE_gauge`,
`phiI_gauge` against the exact field shifted by the volume mean of the `phiE`
error; the gauge rows are the convergence quantities.

### Stack

- myocardium solver: `bidomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov` or `sbdf2` (`-coupling`), `ddt(Vm)` `backward`
  or `Euler` (`-ddt`), `nOuterCorrectors 2`
- ionic model: `bidomainFDAManufactured`
- field verifier: `manufacturedFDABidomainVerifier`
- gradient: `leastSquares` or `Gauss linear` (`-gradient`)
- linear solvers: PCG/DIC, absolute tolerance `1e-15` (`-phiTolerance` for
  `phiE`)

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
| `-dim 1\|2\|3` | dimension of the hex box | 3 |
| `-N n` | cells per unit length; tet: element size `1/n` | 20 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-gradient leastSquares\|GaussLinear` | cell-gradient reconstruction | leastSquares |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |
| `-coupling godunov\|sbdf2` | ionic-PDE time coupling | godunov |
| `-ddt backward\|Euler` | time scheme of `Vm` | backward |
| `-phiTolerance x` | absolute tolerance of the `phiE` solve | 1e-15 |

`postProcessing/<d>D_<n>_cells.dat` holds the field errors (`Vm`, `phiE`,
`phiE_gauge`, `phiI`, `phiI_gauge`, `u1`, `u2` and the bulk, boundary and
gradient rows of
[Bulk, boundary and gradient rows](../README.md#bulk-boundary-and-gradient-rows));
on tetrahedra `n = round(cbrt(cells))`.

Regression (the defaults: 3D hexahedra, `N = 20`, serial): the field errors
and their bulk, boundary and gradient rows.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables.

| study | question |
|---|---|
| [cartesianConvergence](setup/studies/cartesianConvergence/README.md) | spatial order on hexahedra, 1D/2D/3D |
| [tetConvergence](setup/studies/tetConvergence/README.md) | spatial order on tetrahedra, `leastSquares` against `Gauss linear` |
| [temporalConvergence](setup/studies/temporalConvergence/README.md) | temporal order of the four coupling and time-scheme combinations |
| [linearToleranceControl](setup/studies/linearToleranceControl/README.md) | effect of a loose `phiE` tolerance on the tetrahedral errors |
| [tetTemporalControl](setup/studies/tetTemporalControl/README.md) | share of the time step in the tetrahedral errors |

## Variants

- [sphericalShell](sphericalShell/README.md): sealed curved walls crossed by a
  conductivity that rotates through the wall.
