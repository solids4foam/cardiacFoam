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
