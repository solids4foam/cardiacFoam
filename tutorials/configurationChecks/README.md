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

[`comparisons`](comparisons) then checks that an option changes the results
when it should (for example explicit versus implicit monodomain), or leaves
them unchanged when it should (for example adding an ECG, or a field
conductivity equal to the uniform one).

## Layout

| Path | Contents |
| --- | --- |
| `base/slab` | 10 mm x 10 mm 2-D slab (20 x 20 cells) with a small Purkinje graph; monodomain, bidomain and eikonal |
| `base/bath` | 1-D bath-myocardium-bath column (60 cells); bidomain with the extracellular potential domain |
| `base/singleCell` | one cell; `singleCellSolver` |
| `variants/<name>` | files copied over a base, at least `constant/electroProperties` |
| `regression/manifests` | files written by each `pass` check |

Invalid configurations are mostly a valid variant plus a one-entry
`foamDictionary` edit, given in the `checks` file.

## Running

```bash
./Allrun                                  # all checks, about 1-2 minutes
./Allrun bidomain torsoECGWithoutBath     # selected checks
./regression/regressionTest.sh --update-manifests   # after an intended change
./Allclean
```

Each check runs in `runs/<name>`, where its logs and outputs can be
inspected. Serial runs are deterministic, so `comparisons` uses exact file
equality.
