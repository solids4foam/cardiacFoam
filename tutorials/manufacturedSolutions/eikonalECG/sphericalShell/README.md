# eikonalECG/sphericalShell

Eikonal manufactured solution on a spherical heart shell, meshed as a cubed
sphere (`-mesh hex`, default) or with gmsh tetrahedra (`-mesh tet`): sealed
curved boundaries approximated by planar faces, a rotated anisotropic
conductivity with all six tensor components that turns through the wall, and
an activation time whose conormal flux vanishes on both spheres and whose
wavefronts collide on the walls. The case checks the sealed-wall operators and
the outer iteration of the eikonal solver on a field with non-zero normal and
tangential gradients on the walls; with `-ecg` it also checks the eikonal
ECG of that field. It does not check a stimulus region (the verifier pins the
cells near the stimulus point) or a time dependence of the activation time.

Origin: this work.

## Problem

```text
-div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1 + S    in r1 < |x| < r2
n . M grad psi = 0                                            on |x| = r1, |x| = r2
```

with `M = G/(chi Cm)`. The shell has `r1 = 0.5` and `r2 = 1`; the inner
sphere is the patch `endocardium` and the outer sphere the patch
`epicardium`. `G` is the conductivity of the
[monodomain shell](../../monodomainPseudoECG/sphericalShell/README.md), and
`chi = 3`, `Cm = 2`, `c0 = 1`.

`psi` is the smoothed anisotropic distance from a mid-wall stimulus point
`xs = 0.75 (1, 0.3, 0.2)/|(1, 0.3, 0.2)|`,
`sqrt((x - xs)^T M(xs)^-1 (x - xs) + 0.05^2)/c0`, plus the sealing correction
of the monodomain shell; `S` is the exact residual of the operator,
including the gradient of `M`. The sealed problem fixes `psi` up to a
constant, so the verifier pins the cells whose centre lies within `0.15` of
`xs` to the exact `psi` (the nearest cell when the mesh is too coarse to have
a centre in that ball) and imposes no value on a patch. `psi` has critical
points on the walls, where the metric `sqrt(grad psi . M grad psi)` is not
differentiable; the generator lists them
(`python3 shell_eikonal_mms.py --report`, in
[the generator folder](../../../generators/README.md)).

### Stack

- myocardium solver: `eikonalSolver` with `useGraphPrePopulation true` and
  `eikonalAdvectionDiffusionApproach true`: the metric term is the convection
  term `u . grad psi`, `u = c0 M grad psi/|grad psi|_M`, discretised by
  `div(phiU,activationTime) Gauss linearUpwind limitedGrad` with
  `limitedGrad cellMDLimited leastSquares 1` (`system/fvSchemes`);
  `-formulation explicit` evaluates the metric from the cell gradient instead
- verifier: `manufacturedShellEikonalVerifier`
- ECG (`-ecg` only): `eikonalECG` with `manufacturedEikonalECGVerifier`,
  `domain shell`
- conductivity: `conductivitySource field`; `0/Conductivity` written by
  `setManufacturedFields`
- walls: `sealedHeartBoundary true`
- gradient: `sealedWallLeastSquares` for `grad(activationTime)` and
  `grad(activationTime.component(0))` (`system/fvSchemes`): wall cells take no
  value from the sealed walls; `leastSquares` elsewhere
- stimulus box: `stimulusLocationMin/Max` lie outside the shell. With the
  verifier active `applyConstraints` rewrites the activation time of every
  cell before the constraints are collected, so a box inside the shell has no
  effect
- solver controls (`system/fvSolution`): those of
  [eikonalECG](../README.md) (`nOuterCorrectors 8000`, outer residual
  tolerance `1e-8`, equation relaxation `0.9` on every outer iteration, `PBiCGStab`/`DILU` for the
  advection-diffusion matrix), with `PCG`/`DIC` at tolerance `1e-11` for the
  symmetric matrix of the explicit metric; see
  [Outer iteration](#outer-iteration)

## Mesh

- `-mesh hex` (`system/blockMeshDict`): six projected blocks of a cubed
  sphere. `N`, a multiple of 4, sets `k = N/4`: `3k` cells along a cube-face
  edge and `2k` radially, so the shell has `N` radial cells per unit length.
- `-mesh tet` (`shell.geo.template`): Delaunay tetrahedra of size `1/N` from
  `gmsh -setnumber lc 1/N -setnumber seed <seed>`, with one thread, so a gmsh
  version gives the same mesh for the same `N` and seed; imported with
  `gmshToFoam`.

The mesh has exactly two patches, `endocardium` and `epicardium`, and one
cellZone, `myocardium`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/`, `0/Conductivity` and the mesh. `Allrun` stops
with `run ./Allclean first` if a mesh log exists, and, for `-mesh tet`, when
`gmsh` or `gmshToFoam` is missing or gmsh reports an error.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (4 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-mesh hex\|tet` | cubed sphere or gmsh tetrahedra | hex |
| `-N n` | hex: a multiple of 4, radial cells per unit length; tet: element size `1/n` | 12 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-wallGradient sealed\|stock` | `sealedWallLeastSquares` or the stock `leastSquares` for `grad(activationTime)` | sealed |
| `-formulation explicit\|advectionDiffusion` | `eikonalAdvectionDiffusionApproach` false or true | advectionDiffusion |
| `-iterations n` | at most `n` outer iterations | 8000, stopped by the residual control |
| `-exactStart` | start the iteration from the exact field outside the pinned cells (`initialiseExact`); with `-iterations 1` the field after one linearised solve from the exact one | off |
| `-pinnedRadius x` | radius of the ball of pinned cells | 0.15, the header's `rPin` |
| `-ecg` | the eikonal ECG at E1-E5 and its shell reference (writes `local/ecg`) | off |
| `-template singleEvent\|sine` | with `-ecg`: the manufactured ECG template and its sampling (writes `local/ecgTemplate` for `sine`) | singleEvent |

`postProcessing/3D_<hex|tet>_<N>_cells.dat` holds the error summary
(`field L1 L2 Linf`), with the header lines `# h`, `# cellsPerDirection`,
`# pinnedCells` and `# pinnedRadius` and the `activationTime`, bulk, boundary
and gradient rows of
[Bulk, boundary and gradient rows](../../README.md#bulk-boundary-and-gradient-rows).
The outer residual history is the `Initial residual` of the
`activationTime` solves in `log.cardiacFoam`.

## ECG (`-ecg`)

The `ecgDomains` block of `constant/electroProperties` is read only when
`local/ecg` sets `_ecg yes;`. `eikonalECG` reconstructs
`Vm = U(t - psi)` from the computed activation time and sums
`G grad Vm . (x - e)/|x - e|^3 V` over the cells, with
`grad Vm = -U'(t - psi) grad psi` and `grad psi` the cell gradient of
`fvSchemes` (so `-wallGradient` reaches the ECG). The verifier's reference at
electrode `e` and sample time `t` is the same functional of the exact field,

```text
int_{r1 < |x| < r2} -U'(t - psi) (G grad psi) . (x - e)/|x - e|^3 dV
```

with `psi`, `grad psi`, `G`, `r1` and `r2` from
`manufacturedShellEikonalReference.H` and `c0` from `eikonalSolverCoeffs`
(no `sigmaExtracellular`, as in the cube case). The integrand does not
separate in space and time, so the verifier evaluates `psi` and the kernel
`-w (G grad psi) . (x - e)/|x - e|^3` once at the nodes of each order and
forms every sample as the sum of the kernel times `U'(t - psi)`; in parallel
the radial nodes are shared out among the ranks.

Order `q` is a spherical tensor rule: `q` Gauss-Legendre nodes in `r`
(Jacobian `r^2`), `q` in `mu = cos(theta)` and `2q` trapezoid nodes in `phi`.
The smoothed distance of `psi` has a near-conical tip at the stimulus point
(width `0.05/sqrt(lambda)`, 0.002 to 0.008, with `lambda` the eigenvalues of
`M(xs)^-1`), so the rule converges algebraically rather than spectrally:

| template | reference | checks | error of the reference (relative to the trace maximum) |
|---|---|---|---|
| singleEvent | `q = 128` | `24 48 96` | at most `1.3e-5`, at `t = 0` when the front is at the stimulus point |
| sine | `q = 192` | `48 96 128` | at most `6e-4`: `U'` repeats about 28 times across the activation range |

The errors were measured against an independent evaluation (a smooth cut-off
of radius `0.15` around the stimulus point, the inner part in coordinates
centred there with the radius mapped so the smoothed distance is
`0.05 cosh v`, the outer part with this rule at `q = 384`). The singleEvent
reference costs about 40 s and 0.3 GB, the sine reference about 30 s and
1.2 GB, in serial; the `deltaQuadrature` columns compare the checks with the
reference. The verifier stops if an electrode lies on or inside `|x| = 1`,
or if `chi`, `Cm` or the ECG conductivity are not those of the header.

The template, its parameters and the sampling are those of the
[cube case](../README.md#single-event-template): `templateDuration 64`,
`templateRise 4`, sampled from `-6` to `92` by `0.25` (singleEvent), or
`sin(2 pi s)` from `0` to `2` by `0.01` (sine). On the shell `psi` runs from
0.05 to 28.07, so the single event no longer ends inside the window
everywhere: the last-activated points are still repolarising at `t = 92`.

| electrode | position | `\|e\|` |
|---|---|---|
| E1 | (0.4, -0.8, 0.8) | 1.2 |
| E2 | (0.6, 1.2, 0.4) | 1.4 |
| E3 | (-0.96, 0, -1.28) | 1.6 |
| E4 | (-0.6, 1.2, -1.2) | 1.8 |
| E5 | (1.2, -0.96, -1.28) | 2.0 |

Outputs: `postProcessing/eikonalECG.dat`,
`postProcessing/manufacturedEikonalECG_ECG.dat` (per sample: numeric,
`refQ<q>`, `errQ<q>`, `deltaQuadratureQ<q>_Q<qRef>`) and
`postProcessing/manufacturedEikonalECGSummary_ECG.dat` (L1, L2, Linf over the
samples). The field outputs are the same with and without `-ecg`.

```bash
./Allrun -mesh tet -N 10 -wallGradient stock -ecg -template singleEvent
```

## Outer iteration

The eikonal equation is steady. `cardiacFoam` solves it with an outer fixed
point: each outer iteration recomputes the velocity
`u = c0 M grad psi/|grad psi|_M` from the current activation time, assembles
the linear advection-diffusion equation with that velocity, and solves it.
The run stops when the initial residual of an `activationTime` solve falls
below the `residualControl` tolerance `1e-8`, or after `nOuterCorrectors`
iterations. `log.cardiacFoam` ends with `PIMPLE: converged in <n>
iterations` in the first case and `PIMPLE: not converged within <n>
iterations` in the second; the error summary is written in both cases, so
check the log before using a result.

The case relaxes the `activationTime` equation by `0.9` on every outer
iteration (`relaxationFactors/equations` in `system/fvSolution`; the pattern
`"activationTime.*"` also matches the final iteration `activationTimeFinal`).
Equation relaxation adds `(1/0.9 - 1)` times the matrix diagonal to both
sides of the equation, with the current activation time on the right. It
does not change a converged solution, because the added terms cancel when
the activation time stops changing, but it shortens every outer step.

- Why the case relaxes. Without relaxation the outer iteration converges
  quickly on coarse meshes, but on finer meshes it does not settle: the
  residual stops falling and oscillates, or slowly rises, while the
  activation time keeps changing. With relaxation `0.9` every level of the
  studies converges to `1e-8`. Where both converge, the two give the same
  errors to within what the `1e-8` tolerance allows.
- What it costs. Each relaxed step moves the activation time only part of
  the way, so the number of outer iterations grows roughly with the square of
  `N`. The default cap `nOuterCorrectors 8000` is reached on the finer levels
  of the studies; raise it with `-iterations n`, which writes `local/outer`.
  A larger cap only removes the stop; it does not change the solution, and a
  run that converges stops at the tolerance whatever the cap.
- Changing the relaxation. There is no `Allrun` option for it. In a copy of
  the case, edit the factor in `system/fvSolution` with a text editor, or
  delete the `relaxationFactors` block to run without relaxation. Do not use
  `foamDictionary -set` or `-remove` on `system/fvSolution` or
  `constant/electroProperties`: it rewrites the file without its
  `#includeIfPresent` lines, so `-iterations` and `-formulation` no longer
  reach the case. Then watch the residual in `log.cardiacFoam`: a residual
  that stops decreasing well above `1e-8` and keeps fluctuating, or rises, is
  the outer iteration failing to converge, and more iterations will not fix
  it.
- Larger `c0`. The case has `c0 = 1`. The cell Peclet number
  `c0 h/sqrt(m)`, with `m` an eigenvalue of `M`, is at most about 4 on the
  coarsest study meshes and falls below 1 in every direction from tet
  `N = 20` and hex `N = 40`, so diffusion carries much of the equation.
  `c0` is read from `eikonalSolverCoeffs` in `constant/electroProperties`
  (edit it by hand, as above), and the exact solution and its source follow
  it (`psi` scales as `1/c0`), so the same case runs with advection
  dominating. With larger `c0` the outer
  iteration does not reach `1e-8` on the tested meshes, with or without
  relaxation.

The iterations each level needs, the levels where the unrelaxed iteration
fails, and the larger-`c0` runs are in the study READMEs
([hexConvergence](setup/studies/hexConvergence/README.md#outer-iterations-and-relaxation),
[tetConvergence](setup/studies/tetConvergence/README.md#outer-iterations-and-relaxation)).

Regression (the case defaults: hex cubed sphere, `N = 12`, parallel on four
ranks; no gmsh needed): the `activationTime` errors and the bulk, boundary and
wall-gradient rows.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables, with the outer-iteration counts and residuals.

| study | question |
|---|---|
| [hexConvergence](setup/studies/hexConvergence/README.md) | spatial order and outer iteration on the cubed sphere, sealed-wall against stock gradient, one iteration from the exact field, explicit formulation |
| [tetConvergence](setup/studies/tetConvergence/README.md) | spatial order and outer iteration on gmsh tetrahedra, sealed-wall against stock gradient |
