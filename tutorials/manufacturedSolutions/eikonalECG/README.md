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

The activation-time field (`-conductivity rotated`, default, or `axis`; 1D/2D/3D)
and the eikonal-driven ECG are defined in the construction document
[eikonal_mms.pdf](../generators/docs/eikonal/eikonal_mms.pdf), sections
"Cartesian eikonal field" and "Eikonal-driven ECG", with the generators
[`cartesian_eikonal_mms.py`](../generators/cartesian_eikonal_mms.py) and
[`eikonal_ecg_mms.py`](../generators/eikonal_ecg_mms.py) and the NumPy
references [`cartesian_eikonal.py`](../generators/reference/cartesian_eikonal.py)
and [`eikonal_ecg.py`](../generators/reference/eikonal_ecg.py). The exact
activation time is imposed on every face and on the cells of smallest centre
`x`, `y` or `z`. The solver evaluates the ECG with the reconstructed `grad psi`
of the computed field; the reference is a Gauss-Legendre order-48 quadrature of
the exact integrand.

### Single-event template

`template singleEvent;` in the ECG `verificationModel` (the default) is one
action-potential-like event. `templateDuration` and `templateRise`, defaults
`64` and `4`, are set from the activation range `e^3 - 1 = 19.1`:

- `templateDuration` (about the range over 0.3) exceeds the range, so upstroke and
  repolarisation never meet in one sample and `U'(t - psi)` keeps one sign
  across the cube;
- `templateRise` (about 0.2 times the range) keeps the upstroke near one cell wide at
  `x = y = z = 1`, where `|grad psi|` is twenty times its value at the origin,
  on the coarsest mesh.

The sampling window, `-6` to `92` by `0.25`, covers the whole event at every
point.

### Sine template

`template sine;` (`-template sine`) selects the periodic template, sampled from
`0` to `2` by `0.01`. The sine repeats
about 19 times across the activation range, so its ECG integrand changes sign
many times over the cube and the integral cancels. The integrand is also
spatially under-resolved on coarse meshes, so its quadrature error falls
faster than the discretisation order there, and the ECG can appear to converge
while the activation field does not. The studies
([cartesianConvergence](setup/studies/cartesianConvergence/README.md),
[tetConvergence](setup/studies/tetConvergence/README.md)) measure both
templates.

### Stack

- myocardium solver: `eikonalSolver`, the metric term `c0 |grad psi|_M` as
  the bounded convection term `u . grad psi` (`div(phiU,activationTime) Gauss
  linearUpwind limitedGrad`) or evaluated explicitly (`-formulation`)
- outer fixed point: `nOuterCorrectors 8000`, relaxation 0.9, residual
  tolerance `1e-8` (`-tolerance`)
- field verifier: `manufacturedEikonalVerifier`
- ECG: `eikonalECG` with `manufacturedEikonalECGVerifier` (`domain cube`,
  the default, where `sphericalShell/` uses `domain shell`)
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
| `-template singleEvent\|sine` | manufactured ECG template and its sampling | singleEvent |
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
rows, the ECG errors of E1-E5 and their trace values at `t = 10`, on the
upstroke, and `t = 75`, on repolarisation.

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
  rotated conductivity, with a pinned stimulus ball; `-ecg` adds the ECG of
  the shell field.

## Notes on dictionary entries

**Boundary condition of `0/activationTime`.** `0/activationTime` holds `-1`
for "not yet activated". On the cube, the values on every non-empty patch face
are overwritten each step by `manufacturedEikonalVerifier::applyConstraints`
with the exact manufactured solution. Only the patch *type* matters: a
Dirichlet type keeps those values, whereas `zeroGradient` would extrapolate
them away from the interior and silently invalidate the verification. On the
spherical shell both spheres are sealed walls: the verifier pins the cells
near the stimulus point and imposes no value on a patch.

**Electrodes.** The cube uses the 161 electrodes of `monodomainPseudoECG`:
`E1-E5`, then `R1-R156` with a direction uniform on the sphere and a standoff
uniform in `[0.05, 1.0]`. The spherical shell uses five directions with
`|e| = 1.2, 1.4, 1.6, 1.8, 2.0`, as in `monodomainPseudoECG/sphericalShell`.

**Spherical shell.** The mesh is the cubed-sphere shell described in
`../bidomain/README.md`. The eikonal-driven ECG of the shell field is off
unless `Allrun -ecg` writes `ecg yes;` to `local/ecg`; the template, its
parameters and its sampling are those of the cube case, and `-template sine`
writes `local/ecgTemplate`. The reference integrates the same kernel over the
shell with a spherical Gauss rule; electrodes lie outside the outer sphere
`|x| = 1`.
