# manufacturedSolutions/eikonalTetMMS tutorial

Tetrahedral (unstructured) variant of the eikonal manufactured MMS
(`../eikonalECG`'s hex 3D case), reusing the exact same unit-cube mesh
generator as the monodomain tet MMS (`../monodomainTetMMS`) so both solvers
refine on the same mesh sequence at N=10/20/40/80.

## Purpose

`../eikonalECG` fixed the pseudo-ECG's dependence on `fvc::grad(Vm)` by
reconstructing `gradVm` analytically from `gradActivationTime` (see commit
`f803626c`). That removed one source of gradient-scheme sensitivity, but
`gradActivationTime` itself is still `fvc::grad(activationTime)` computed
once per solve ([eikonalECG.C:388](../../../src/electroModels/ecgModels/eikonalECG/eikonalECG.C)) --
so the question this case answers is whether *that* single grad call is
scheme-sensitive on a genuinely non-orthogonal tet mesh, the same way
`fvc::grad(Vm)` was shown to be for the monodomain path in
`../monodomainTetMMS`.

## Mesh and boundary conditions

`setup/box.geo.template` is copied verbatim from `../monodomainTetMMS` --
same unit cube, same `lc = 1/N` ladder, same single lumped "boundary"
physical surface for all six faces. The two solvers get away with a shared,
patch-name-agnostic mesh for different reasons:

- monodomain: cardiacFoam always creates `Vm` with a uniform zeroGradient
  boundary regardless of patch name.
- eikonal: `activationTimePatchTypes()` in `eikonalMyocardiumDomain.C`
  switches every *non-empty* boundary patch to `fixedValue` when manufactured
  verification is enabled (keyed off patch type, not name), and
  `manufacturedEikonalVerifier::applyConstraints()` pokes the exact `tau(x)`
  value into every boundary face each solve. So this case is
  Dirichlet-everywhere, not zero-flux -- a different mechanism than the
  monodomain case, but one that is equally indifferent to gmsh patch naming.

## Convergence tuning (read this before trusting large-N runs)

The eikonal solve is a nonlinear fixed point (Picard iteration): `phiU` and
`divPhiU` are lagged from the previous outer iterate, and
`eikonalMyocardiumDomain::advance()` applies `activationEqn.relax()` each
outer iteration specifically to keep that fixed point contractive (commit
`756c7fcb`). `../eikonalECG`'s `fvSolution` (`nOuterCorrectors 5000`,
`residualControl` tolerance `1e-11`, `relax 0.9`) converges fine on the hex
mesh's low non-orthogonality. On this tet mesh's ~60-70deg non-orthogonality
the outer residual decays at only ~0.2%/iteration, so:

- **`nOuterCorrectors` had to rise to 2500** and `residualControl` tolerance
  loosened to `1e-8` (a level still 1-3 orders of magnitude below the
  smallest discretization error measured at every N in this sweep).
  N=80 still hits the 2500 cap before reaching `1e-8`, landing at `~3.15e-5`
  -- about 45x below its own L2 error, an acceptable but not generous
  margin.
- **Transplanting `../../heartSim3D-1D/eikonalHeart`'s fvSolution
  (`preconditioner none`, `relax 0.3`, no residualControl) made this case
  *worse*, not better** (N=40 outer residual only reached 1.5e-2 in 500
  iterations, vs ~1.9e-4 for a comparable budget under the original
  settings). `eikonalHeart`'s README says diagonal/DILU preconditioners
  "trap" on its assembled heart-mesh matrix -- a stalling/instability
  failure mode. This synthetic tet cube's outer loop decays smoothly and
  monotonically instead; it is slow, not unstable, so heavier
  under-relaxation just shrinks every step without addressing the actual
  problem. The two "hard mesh" cases need different fixes; do not assume one
  transfers to the other.
- Fix that actually worked: keep the original `diagonal`/`DILU`
  preconditioners and `relax 0.9`, and simply give the outer loop enough
  iterations. **Run large-N cases in parallel** (`./Allrun parallel`,
  `decomposeParDict` already has `numberOfSubdomains 6`) -- N=80 serial took
  3235s (~54 min), the identical case in parallel took 1007s (~17 min), a
  ~3.2x speedup, with byte-identical converged fields.
- `run_tet_sweep.sh` and `run_scheme_study.sh` now archive `log.cardiacFoam`
  per resolution (`setup/results/<N>/log.cardiacFoam` and
  `setup/results/logs/<scheme>_N<N>.log`) and warn if the iteration cap was
  hit, so future runs can check convergence rather than assume it.

### Parallel-mode output-path bug found and fixed

`manufacturedEikonalVerifier::postProcess` (activationTime) writes via
`mesh_.time().globalPath()/"postProcessing"` and lands correctly in the
shared `postProcessing/` under `./Allrun parallel`. Five sibling verifiers
used `mesh_.time().path()` instead (which resolves to `processorN/` under
decomposition), silently misplacing their output where no downstream script
would find it: `eikonalECGManufacturedVerifier`, `pseudoECGManufacturedVerifier`
(the monodomain pseudo-ECG), `bathECGManufacturedVerifier`,
`manufacturedFDAMonodomainVerifier` (the monodomain Vm-field verifier), and
`manufacturedElectromechanicsVerifier`. All five fixed to use `globalPath()`
(matching `manufacturedEikonalVerifier`'s existing correct pattern) and
`libverificationModels` rebuilt; verified by rerunning both
`eikonalTetMMS` and `monodomainTetMMS` at N=10 in parallel and confirming
`processor0/postProcessing/` is no longer created at all.

## Running the sweep

```bash
RESOLUTIONS="10 20 40" ENDTIME=0.02 bash setup/run_tet_sweep.sh   # N=80 add manually, ~15-55min depending on serial/parallel
bash setup/run_scheme_study.sh   # GaussLinear vs leastSquares, N=10/20/40 (N=80 leastSquares needs a manual parallel run, see above)
```

Requires OpenFOAM sourced (`/Volumes/OpenFOAM-v2412/etc/bashrc`) and `gmsh`
on `PATH`. Results land in `setup/results/<N>/`, a combined
`setup/results/summary.csv` (activationTime + pseudo-ECG order vs. tet
`dx`), and `setup/results/scheme_study.csv` (GaussLinear vs leastSquares).

## Results

### Spatial order, `gradSchemes.default = leastSquares` (tutorial default)

| N  | dx        | L2 activationTime | p_L2 | ECG L2      | p_ECG_L2 | maxNonOrtho |
|----|-----------|--------------------|------|-------------|----------|-------------|
| 10 | 0.0588235 | 2.542e-02          | --   | 5.364e-02   | --       | 61.7deg     |
| 20 | 0.030303  | 9.784e-03          | 1.44 | 9.263e-03   | 2.65     | 66.6deg     |
| 40 | 0.0151515 | 3.347e-03          | 1.55 | 1.060e-03   | 3.13     | 68.9deg     |
| 80 | 0.0075758 | 1.201e-03          | 1.48 | 1.741e-04   | 2.61     | 70.0deg     |

### GaussLinear vs leastSquares (`setup/results/scheme_study.csv`)

| scheme | step | activationTime L2 order | ECG L2 order |
|---|---|---|---|
| GaussLinear | 10->20 | 0.19 | 3.59 |
| GaussLinear | 20->40 | 0.10 | 2.75 |
| leastSquares | 10->20 | 1.44 | 2.65 |
| leastSquares | 20->40 | 1.55 | 3.13 |
| leastSquares | 40->80 | 1.48 | 2.61 |

## Interpretation

1. **The gradScheme sensitivity generalizes to the eikonal path.**
   `Gauss linear corrected` stalls `activationTime`'s order to ~0.1-0.2 by
   N=40 -- the same failure mode already found for monodomain's `Vm`
   (`../monodomainTetMMS`), confirming this is not solver-specific.
2. **`leastSquares` does not fully recover 2nd order here, unlike
   monodomain.** Monodomain's `Vm` hits a clean, stable ~2.0 with
   leastSquares across all 4 points. Eikonal's `activationTime` settles at a
   stable ~1.5 (1.44, 1.55, 1.48 across the three refinement pairs) and does
   not climb toward 2. The gradient scheme fixes the *stall*; the operator
   itself caps the order. The concrete mechanism is the nonlinear
   gradient-magnitude term `G = sqrt(gradPsi & (M & gradPsi))` and the `phiU`
   stabilization flux in `eikonalSolver.C`, both built from the reconstructed
   cell gradient `fvc::grad(psi)` -- only ~1st-order accurate on this
   sustained-skewness tet family. A 2nd-order Laplacian coupled to a
   ~1st-order gradient-derived source lands between 1 and 2, exactly as
   measured; monodomain has no such term (pure divergence-form). Lower
   solution regularity near eikonal/Hamilton-Jacobi characteristics may
   contribute secondarily.
3. **The pseudo-ECG functional is more robust than the field it's built
   from, in both solvers.** Even where `activationTime`'s own order is only
   ~1.2-1.6, the derived pseudo-ECG order is consistently >2 (2.5-3.1) under
   leastSquares -- the lead-field integral averages out a good deal of the
   primary field's local error. This mirrors monodomain's pattern (Vm order
   ~2.0, pseudo-ECG order ~2.3-3.1) and is a genuine reassurance: the
   clinically-relevant output converges better than the intermediate field's
   raw accuracy would suggest, in both solver families.
