# Outer-corrector study for the staggered monodomain solve

## Question

Does repeating the implicit monodomain diffusion solve through
`nOuterCorrectors` materially change the solution on a non-orthogonal
tetrahedral mesh, and is the configured `nNonOrthogonalCorrectors` entry active
in the current code path?

This study deliberately does **not** restore the removed
`pimpleStaggeredElectrophysicsAdvanceScheme`. Commit `02d5ca87` removed that
scheme because nested ownership of `pimpleControl::loop()`, repeated ODE time
advances, and accumulating PVJ sources made it an invalid strong-coupling
algorithm. Every run here uses the remaining single-pass
`staggeredElectrophysicsAdvanceScheme`. Only the implicit tissue diffusion
solve is repeated.

## Hypotheses

1. `nOuterCorrectors` can affect a corrected Laplacian because each outer sweep
   reassembles the explicit non-orthogonal/cross-diffusion contribution using
   the newest `Vm` field.
2. The effect should converge as the outer count increases. The smallest
   acceptable count is the first count whose field difference from the
   highest-count reference is at least two orders of magnitude below the MMS
   discretisation error.
3. `nNonOrthogonalCorrectors` now invokes the equation-level
   `correctNonOrthogonal()` loop. At fixed `nOuterCorrectors`, a value $k$
   should produce $k+1$ `Vm` solves per outer sweep and should update the
   deferred cross-diffusion source from the latest field.

## Experiments

The runner performs two experiments on the same generated tetrahedral mesh and
time controls:

| Experiment | Varied setting | Fixed setting | Purpose |
|---|---|---|---|
| `outer` | `nOuterCorrectors = 1,2,3,4,8` | `nNonOrthogonalCorrectors = 0` | convergence and cost of outer reassembly sweeps |
| `nonorth` | `nNonOrthogonalCorrectors = 0,1,2` | `nOuterCorrectors = 1` | verify solve counts and measure convergence of the deferred correction |

For every run it records:

- MMS `Vm` L2 and Linf errors;
- L2 and Linf field differences relative to the experiment reference;
- number of `Vm` linear solves and physical timesteps;
- execution and wall-clock time reported by OpenFOAM;
- SHA-256 of the final ASCII `Vm` field;
- the exact `fvSolution`, `controlDict`, log, and verifier summary.

The high-count outer run is a numerical reference for *corrector convergence*,
not a reference solution of the PDE.

## Running

Quick smoke study (four physical steps on the `N=10` mesh):

```bash
bash setup/correctorStudy/run_corrector_study.sh
```

Paper study with three timing repeats and the spatial ladder:

```bash
RESOLUTIONS="10 20 40" \
ENDTIME=0.2 \
REPEATS=3 \
bash setup/correctorStudy/run_corrector_study.sh
```

Useful overrides:

- `OUTER_COUNTS="1 2 3 4 8"`
- `NONORTH_COUNTS="0 1 2"`
- `NONORTH_FIXED_OUTER=1`
- `OPENFOAM_BASHRC=/Volumes/OpenFOAM-v2412/etc/bashrc`
- `RESULTS_DIR=/absolute/output/path`
- `KEEP_WORK=1`

By default, output is written to `setup/correctorStudy/results/` and temporary
case state to `/tmp/cardiacfoam-corrector-study-*`. The source tutorial is not
cleaned or modified.

## Interpretation and decision rule

For each resolution, let `delta_k` be the L2 difference between the field from
`k` outer sweeps and the highest-count reference, and let `E_MMS` be the L2 MMS
error of that reference.

- Accept `k` sweeps as converged for the tested case when
  `delta_k / E_MMS <= 0.01` and the next sweep count changes the MMS error by
  no more than 1%.
- If `k=1` already satisfies the rule at every resolution, remove the outer
  repetition from the monodomain path (or set its documented default to one).
- If more sweeps are needed, report the convergence evidence and retain the
  smallest passing count; do not describe it as a strong 1D--3D coupling loop.
- The `nonorth` experiment must produce
  `nSteps * nOuterCorrectors * (nNonOrthogonalCorrectors + 1)` linear solves.
  Select a non-orthogonal count from convergence of the field change, not from
  the dictionary value alone.

The bath-bidomain PDE split is studied separately under
`bathBidomainTetMMS/setup/couplingStudy`: the global `phiE` correction and the
`phiE`/`Vm` coupling are distinct controls, and a count selected for the
myocardium equation must not be generalized to them without evidence.

## Production decision

Post-fix comparisons at `N=10,20` showed that
`nOuterCorrectors=2, nNonOrthogonalCorrectors=0` and
`nOuterCorrectors=1, nNonOrthogonalCorrectors=1` produce bitwise-identical
fields. The reported manufactured cases use the latter configuration: one
outer sweep owns the physical equation sequence, while one additional
non-orthogonal assembly owns the deferred cross-diffusion update.
