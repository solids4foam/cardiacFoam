# monodomain/cube

## What it does

Monodomain manufactured solutions on the unit cube, on Cartesian hexahedra
(`-mesh hex`, default) or gmsh Delaunay tetrahedra (`-mesh tet`): the FDA
cosine field with a diagonal conductivity, or an anisotropic sin^2 field with
a full rotated conductivity. The case writes the field errors against the
exact solution and, in 3D, the pseudo-ECG errors against its reference
integral.

- **Exercises:** the interior operators of `monodomainSolver` (diffusion,
  cell-gradient reconstruction, non-orthogonal correction on tetrahedra, ionic
  coupling, time integration) and `pseudoECG`. Both fields have a zero normal
  gradient on every face, so no wall treatment is tested; sealed walls are
  checked by [sphericalShell](../sphericalShell/README.md). Only this case's
  regression checks the pseudo-ECG against an exact reference.
- **Use it for:** verifying a change to the monodomain operators, the time
  coupling or the pseudo-ECG; starting a spatial or temporal convergence study
  (`setup/studies/`); checking a new mesh type against an exact solution
  (`-mesh tet` shows how a gmsh mesh enters the case).

## Regression

Tier: `standard`.

`regression/regressionTest.sh` first checks that `constant/electroProperties`
names the electrodes E1-E5, then runs `./Allrun parallel` with the case
defaults only (3D hexahedra, anisotropic field, `N = 20`, `T = 0.2`,
`godunov` and `backward`, six ranks with scotch); none of the other `Allrun`
options is run. It compares 34 rows: `regression/cube.reference` (19 rows)
holds the cells per direction, the final time, 12 error norms of `Vm`, `u1`,
`u2` and the boundary and gradient rows, and the final pseudo-ECG values of
E1-E5; `regression/cube.ecg.reference` (15 rows) holds the L1, L2 and Linf
pseudo-ECG errors of E1-E5 against the order-96 reference integral. The
expected values are frozen from an earlier run, with absolute tolerances of
about `1e-5` of each value, so the check fails on any change of the results,
not only on a loss of accuracy. No marker applies.

```bash
./regression/regressionTest.sh
```

## Problem

Origin: the FDA field is the monodomain problem of the FDA regulatory science
tool (Pathmanathan & Gray 2014); the anisotropic field and the pseudo-ECG
problem are this work.

The FDA field (`-field fda`, 1D/2D/3D) and the anisotropic field (`-field
anisotropic`, 3D, default) are defined in the construction document
[monodomain_mms.pdf](../../generators/docs/monodomain/monodomain_mms.pdf), sections
"FDA cube", "Rotated anisotropic cube" and "Pseudo-ECG", with the generators
[`fda_monodomain_mms.py`](../../generators/fda_monodomain_mms.py),
[`rotated_monodomain_mms.py`](../../generators/rotated_monodomain_mms.py) and
[`pseudo_ecg_mms.py`](../../generators/pseudo_ecg_mms.py) and the NumPy references
[`fda_monodomain.py`](../../generators/reference/fda_monodomain.py),
[`rotated_monodomain.py`](../../generators/reference/rotated_monodomain.py) and
[`pseudo_ecg.py`](../../generators/reference/pseudo_ecg.py). The pseudo-ECG at 161
electrodes is compared with a Gauss-Legendre order-96 reference integral of the
same field. `sigmaExtracellular` is `1/(4 pi)` = `0.07957747154594767`, for which
`1/(4 pi sigmaExtracellular)` is exactly 1 in double precision, so the ECG is
the unscaled sum the reference integrates.

### Stack

- myocardium solver: `monodomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov` or `sbdf2` (`-coupling`), `ddt(Vm)` `backward`
  or `Euler` (`-ddt`)
- ionic model: `monodomainFDAManufactured`
- field verifier: `manufacturedFDAMonodomainVerifier` or
  `manufacturedAnisotropicMonodomainVerifier` (`-field`)
- pseudo-ECG: `pseudoECG` with `manufacturedPseudoECGVerifier`, its reference
  matched to the field (`anisotropic yes` for the anisotropic field; `domain
  cube`, the default, where `sphericalShell/` uses `domain shell`); 3D only
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
[Bulk, boundary and gradient rows](../../README.md#bulk-boundary-and-gradient-rows));
on tetrahedra `n = round(cbrt(cells))`. In 3D,
`postProcessing/pseudoECG.dat` holds the electrode traces and
`postProcessing/manufacturedPseudoECGSummary_ECG.dat` their errors against
the reference integral.

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

- [sphericalShell](../sphericalShell/README.md): sealed curved walls crossed by a
  conductivity that rotates through the wall.

## Notes on dictionary entries

**Slab meshes.** `system/blockMeshDict.2D` is an `nx x ny x 1` slab, thin in z
with `empty` faces; `system/blockMeshDict.3D` is a unit slab with
`(nx ny nz)` cells in x, y and z. `sphericalShell/system/blockMeshDict` builds
the heart shell `r1 < |x| < r2` as six projected blocks of a cubed sphere;
`N` (a multiple of 4) sets `k = N/4`: `3k` cells along a cube-face edge and
`2k` radially, i.e. `N` cells per unit length radially.

**Electrodes (cube).** The five-point set is kept for the tet per-electrode
table (`tbl-tet-pseudoecg-electrodes`). The random scatter adds 156 points
(`R1-R156`): direction uniform on the sphere, standoff uniform in
`[0.05, 1.0]` from the nearest cube surface point (below 0.05 risks the
near-boundary mesh layer; above ~1.0 the field has decayed too far for the
distance to discriminate). `R1-R56` and `R57-R156` are independent draws,
concatenated so that the first 56 stay fixed if the set is extended. Each
point is placed by exiting the cube along the ray from the centre
`(0.5, 0.5, 0.5)` and going the standoff distance further.

**Electrodes (spherical shell).** Five directions with `|e| = 1.2, 1.4, 1.6,
1.8, 2.0`. The pseudo-ECG of the shell field is off unless `Allrun -ecg`
writes `ecg yes;` to `local/ecg`. The reference integrates
`(G grad F).(x - e)/|x - e|^3` over the shell with a spherical Gauss rule;
electrodes lie outside the outer sphere `|x| = 1`. `Allrun -electrodes <file>`
copies the file in (`name (x y z);` lines, `|x| > 1`).

**Solver settings.** In `system/fvSolution`, one additional deferred-correction
assembly gives the same result as two outer sweeps and assigns the iteration
count correctly.
