# configurationChecks

## What it does

Short runs that check how cardiacFoam turns `constant/electroProperties` into
models: which dictionary configurations are accepted, what each writes, and
which are rejected and with which message. These are not a physical
tutorial. Each line of [`checks`](checks) runs one configuration on a small
base case; [`comparisons`](comparisons) then compares the outputs of pairs of
checks.

- **Exercises:** run-time selection and keyword reading of
  `electroProperties` for the monodomain, bidomain, eikonal and single-cell
  myocardium models, Purkinje networks, the bath and the ECG and lead-field
  options; the set of files each configuration writes; the fatal-error
  messages for unknown selector values, invalid combinations and misspelt
  keywords; parallel against serial runs; restarts from the latest time; and
  the post-processing utilities run on a result. Only this case checks the
  rejection messages and the written-file sets.
- **Use it for:** checking a dictionary change, a new keyword or a new model
  selection before running a real case; after an intended change to what a
  configuration writes, rewriting the manifests with `--update-manifests`.

## Regression

Tier: `standard`.

`regression/regressionTest.sh` (which `Allrun` calls) has its own driver and
message format; it uses only helpers from the shared library, not its
reference grammar. It runs the 135 checks in `checks`, serial except the 18
`np=2` checks, then the 49 comparisons in `comparisons`:

- 64 `pass` checks: every log must end cleanly with no FOAM FATAL error, and
  the written files must match `regression/manifests/<name>` (64 manifests).
  9 of them also run a script from `post/` (`restartLatestTime` 6 times,
  `recomputePseudoECG`, `eikonalECGGradientCheck`,
  `bathBidomainInterfaceMetrics`).
- 71 `fail:<text>` checks: the run or its preparation must stop with a FOAM
  FATAL error that mentions `<text>`.
- 49 comparisons: 13 `same` (identical file), 15 `differ` (the option has an
  effect) and 21 `close:<rtol>:<atol>` (parallel against serial, and the
  recomputed against the solver's pseudo-ECG, with
  `applications/scripts/cardiacCompareCases`; 20 at rtol 1e-6, atol 1e-12,
  one at rtol 1e-7, atol 1e-15). A comparison is skipped unless both of its
  checks ran.

24 checks carry `needs=solids4foam` and are skipped in the `lightweight`
build mode; 27 carry `skip=v2312` and are skipped on OpenFOAM v2312. The case
has no marker files. Each check prints `PASS:`, `FAIL:` or `SKIP:`, and the
run ends with `Configuration checks PASSED (<n> checks, <k> skipped)` or
`Configuration checks FAILED (<f>/<n>): <names>`, where `<n>` counts the
checks and comparisons that ran. The script accepts check names and
`--update-manifests`; it does not accept the shared `--check-only` and
`--report` options. Run it on its own from the case root:

```bash
./regression/regressionTest.sh
```

## What is checked

Each line of [`checks`](checks) is one short run of a configuration on a small
base case:

- `pass` checks must complete cleanly and write exactly the files listed in
  `regression/manifests/<name>`, which catches a field or output file that is
  no longer written, or a new one.
- `fail:<text>` checks must stop with a FOAM FATAL error that mentions
  `<text>`, in `cardiacFoam` or in the case preparation (`Allprepare`). They
  cover unknown selector values, invalid combinations such as `torsoECG`
  without a torso, and misspelt or unknown keywords and values, whose error
  must name the closest valid one (`did you mean ...`).
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
| `base/leadFieldBox` | 3-D heart+torso box (16 x 5 x 5 cells) sheared in y and z, so that every face is non-orthogonal; the heart touches outer patches and its wall is the processor boundary of a 2-subdomain simple decomposition. Its `Allprepare` builds the lead fields; the variants exercise the discrete adjoint with rotated conductivities |
| `base/singleCell` | one cell; `singleCellSolver` |
| `variants/<name>` | files copied over a base, at least `constant/electroProperties` |
| `post/<script>` | post-processing run after the solver by `post=<script>` checks |
| `regression/manifests` | files written by each `pass` check |

Invalid configurations are mostly a valid variant plus a one-entry
`foamDictionary` edit, given in the `checks` file. Edits are applied after
the case preparation, so a variant's own `Allprepare` (for example the
`buildLeadFields` step of the lead-field variants) prepares the unedited
configuration. A variant's `Allprepare.parallel`, if present, runs after
`decomposePar` in `np=N` checks.

## Running

```bash
./Allrun                                  # all checks, about 2 minutes
./Allrun bidomain torsoECGWithoutTorso    # selected checks
./regression/regressionTest.sh --update-manifests   # after an intended change
./Allclean
```

Each check runs in `runs/<name>`, where its logs and outputs can be
inspected. Serial runs are deterministic, so `same` and `differ` comparisons
use exact file equality; parallel runs differ from serial ones at round-off
level, so they are compared within a tolerance.

## Notes on dictionary entries

- `0/activationTime` holds `-1` for "not yet activated"; only the eikonal
  myocardium reads it.
- `0/t` is the transmural coordinate for `ionicHeterogeneity` (0 endocardium,
  1 epicardium).
- `0/Conductivity` is read only when `conductivitySource` is `field`. In
  `monodomainPeriodicFieldConductivity` the fibre conductivity rises from 0.02
  to 0.06 across the 20 columns of cells, so it jumps across the periodic
  faces.
- The `eikonalECG` variant builds its transmural waveform templates from the
  block in its `constant/electroProperties`.
- The single-cell variants (`singleCellBuenoOrovio`,
  `singleCellTWorldActiveTension`) give one beat at 20 ms; `singleCellStimulus`
  times are in ms.
- `base/slab` is a 10 mm x 10 mm x 0.1 mm slab with 20 x 20 x 1 cells
  (0.5 mm); `monodomainPeriodicFieldConductivity` is periodic in x.
  `system/topoSetDict` defines the left half of the slab as a cellZone for the
  cellZone-subset check.
- `base/slab/constant/purkinjeGraph`: root 0 at (1, 5) mm branches at node 1
  (3, 5) mm to the junctions 2 (5, 3) mm and 3 (5, 7) mm, subdivided into
  ~0.5 mm segments so `monodomain1DSolver` has a cable to integrate. The
  junctions are close to the root so that the network activates the slab
  within the short run. `conductionEdges` entries are
  `(nodeA nodeB length conductance)`.
