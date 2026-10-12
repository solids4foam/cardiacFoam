# monodomainPseudoECG

Monodomain manufactured solutions on the unit cube, on Cartesian hexahedra
(`-mesh hex`, default) or gmsh Delaunay tetrahedra (`-mesh tet`), with the
pseudo-ECG checked against its reference integral in 3D. Two fields: the FDA
cosine field with a diagonal conductivity, and an anisotropic sin^2 field with
a full rotated conductivity. Both have a zero normal gradient on every face,
so these cases check the interior operators (diffusion, cell-gradient
reconstruction, non-orthogonal correction, ionic coupling, time integration)
and the pseudo-ECG, not a wall treatment; sealed walls are checked by
[sphericalShell](sphericalShell/README.md).

Origin: the FDA field is the monodomain problem of the FDA regulatory science
tool (Pathmanathan & Gray 2014); the anisotropic field and the pseudo-ECG
problem are this work.

## Problem

```text
chi Cm dVm/dt = div(G grad Vm) - chi Cm Iion + I_s    in [0,1]^d
n . grad Vm = 0                                       on every face
```

`Vm = sqrt(1 + t) F` with

- FDA (`-field fda`, 1D/2D/3D): `F = cos(pi x)`, `cos(pi x) cos(2 pi y)`,
  `cos(pi x) cos(2 pi y) cos(3 pi z)`; `G = diag(0.111453302, 0.121585420,
  0.030396355)`
- anisotropic (`-field anisotropic`, 3D, default): `F = sin^2(pi x)
  sin^2(pi y) sin^2(pi z)`; `G` the full tensor `(0.106875 -0.0084931
  -0.022561 0.116682 -0.0130256 0.0398782)`

with `chi = 3`, `Cm = 2` and the FDA ionic closure (`u1`, `u2`, `u3`). The
pseudo-ECG at 161 electrodes is compared with a Gauss-Legendre order-96
reference integral of the same field.

### Stack

- myocardium solver: `monodomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov` or `sbdf2` (`-coupling`), `ddt(Vm)` `backward`
  or `Euler` (`-ddt`)
- ionic model: `monodomainFDAManufactured`
- field verifier: `manufacturedFDAMonodomainVerifier` or
  `manufacturedAnisotropicMonodomainVerifier` (`-field`)
- pseudo-ECG: `pseudoECG` with `manufacturedPseudoECGVerifier`, its reference
  matched to the field (`anisotropic yes` for the anisotropic field); 3D only
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
| `-dim 1\|2\|3` | dimension of the hex box; the pseudo-ECG runs in 3D only | 3 |
| `-N n` | cells per unit length; tet: element size `1/n` | 20 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-field anisotropic\|fda` | manufactured field, conductivity, verifiers | anisotropic |
| `-gradient leastSquares\|GaussLinear` | cell-gradient reconstruction | leastSquares |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |
| `-coupling godunov\|sbdf2` | ionic-PDE time coupling | godunov |
| `-ddt backward\|Euler` | time scheme of `Vm` | backward |

`postProcessing/<d>D_<n>_cells.dat` holds the field errors (`Vm`, `u1`, `u2`
and the bulk, boundary and gradient rows of
[Bulk, boundary and gradient rows](../README.md#bulk-boundary-and-gradient-rows));
on tetrahedra `n = round(cbrt(cells))`. In 3D,
`postProcessing/pseudoECG.dat` holds the electrode traces and
`postProcessing/manufacturedPseudoECGSummary_ECG.dat` their errors against
the reference integral.

Regression (the defaults: 3D anisotropic field, `N = 20`, parallel on six
ranks): field errors, the final traces of E1-E5 and their pseudo-ECG errors.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables.

| study | question |
|---|---|
| [cartesianConvergence](setup/studies/cartesianConvergence/README.md) | spatial order on hexahedra, 1D/2D/3D, both fields, pseudo-ECG |
| [tetConvergence](setup/studies/tetConvergence/README.md) | spatial order on tetrahedra, `leastSquares` against `Gauss linear` |
| [temporalConvergence](setup/studies/temporalConvergence/README.md) | temporal order of the four coupling and time-scheme combinations |
| [tetTemporalControl](setup/studies/tetTemporalControl/README.md) | share of the time step in the tetrahedral errors |

## Variants

- [sphericalShell](sphericalShell/README.md): sealed curved walls crossed by a
  conductivity that rotates through the wall.
