# configurationChecks

Short runs that check how cardiacFoam turns `constant/electroProperties` into
models. These are not a physical tutorial. They pin down, before the
solids4foam-alignment refactor, which dictionary configurations are accepted,
what each writes, and which are rejected.

## What is checked

Each line of [`checks`](checks) is one short run of a configuration on a small
base case:

- `pass` checks must complete cleanly and write exactly the files listed in
  `regression/manifests/<name>`, which catches a field or output file that is
  no longer written, or a new one.
- `fail:<text>` checks must stop with a FOAM FATAL error that mentions
  `<text>`. They cover unknown selector values and invalid combinations, such
  as `torsoECG` without a bath.
- `np=2` checks run a variant on two processors (simple decomposition) and
  reconstruct it.
- `post=<script>` checks also run a post-processing utility on the result,
  from `post/<script>`, e.g. `recomputePseudoECG`, or restart the run from
  its latest time to twice its end time (`restartLatestTime`).

[`comparisons`](comparisons) then checks that an option changes the results
when it should (for example explicit versus implicit monodomain), or leaves
them unchanged when it should (for example adding an ECG, or a field
conductivity equal to the uniform one). It also checks that every parallel
run matches its serial run up to round-off, and that `recomputePseudoECG`
reproduces the pseudo-ECG written by the solver; these use
`applications/scripts/cardiacCompareCases`.

## Layout

| Path | Contents |
| --- | --- |
| `base/slab` | 10 mm x 10 mm 2-D slab (20 x 20 cells) with a small Purkinje graph; monodomain, bidomain and eikonal |
| `base/bath` | 1-D bath-myocardium-bath column (60 cells); bidomain with the extracellular potential domain |
| `base/singleCell` | one cell; `singleCellSolver` |
| `variants/<name>` | files copied over a base, at least `constant/electroProperties` |
| `post/<script>` | post-processing run after the solver by `post=<script>` checks |
| `regression/manifests` | files written by each `pass` check |

Invalid configurations are mostly a valid variant plus a one-entry
`foamDictionary` edit, given in the `checks` file.

## Running

```bash
./Allrun                                  # all checks, about 2 minutes
./Allrun bidomain torsoECGWithoutBath     # selected checks
./regression/regressionTest.sh --update-manifests   # after an intended change
./Allclean
```

Each check runs in `runs/<name>`, where its logs and outputs can be
inspected. Serial runs are deterministic, so `same` and `differ` comparisons
use exact file equality; parallel runs differ from serial ones at round-off
level, so they are compared within a tolerance.

Lead-field checks use `base/leadFieldColumn`: a bounded passive conductor with
an embedded myocardial cellZone. They validate reciprocal-current scaling and
Wilson weights against a piecewise-linear 1D solution, mono/bidomain projection,
unchanged Vm, canonical/legacy conductivity inputs, MPI output and explicit
heart-only export. The isolated `post/leadFieldGuards` hook also checks stale
inputs, field dimensions, independent sampling and checked restart append.
The column uses a conservative Gauss gradient for its endpoint-voltage oracle.
Separate 3D manufactured tutorials test the vector functional and refinement.

`bidomainLeadFieldRecovery` and its MPI variant prescribe a fixed, nonzero
Vm in a verification-owned heart, then converge the production passive
extracellular recovery and sample with the real `torsoECG` kernel. The same
balanced weights are applied to those electrode potentials. Both paths must
match the independent discrete column voltage .00475 V, and Vm must remain
unchanged. This is a bounded scalar/orthogonal check, not an exact-discrete
adjoint claim for general meshes. Ordinary bathTorsoECG checks exercise the
runtime coupled route separately.
