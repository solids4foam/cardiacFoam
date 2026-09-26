# corrector - bidomain

## Purpose

This study validates the predictor-corrector inner loop convergence and stability for the bidomain solver.

## Execution

```bash
driverFoam sweep-run --spec tutorials/manufacturedSolutions/bidomain/setup/studies/corrector/sweep_corrector_study.json
```

Sweeps `N = 10, 20, 40` across the four reported variants (`baseline`,
`outer2`, `nonorth1`, `combined`), each a `(nOuterCorrectors,
nNonOrthogonalCorrectors)` pair applied as direct
`system/fvSolution:PIMPLE.nOuterCorrectors`/`PIMPLE.nNonOrthogonalCorrectors`
study keys against the `manufacturedBidomain` tutorial record; the short,
fixed step-count screening window (2/9/36 steps) is set via the direct keys
`system/controlDict:writeControl`/`writeInterval`/`writeFormat`.

Comparing the 12 results across variants/resolutions is handled by the
postprocessing module.

**Corrected 2026-09-26 (tutorials-are-pointers, 5.4b-B):** this used to
describe `archive_dir_name`-based output archiving into each case's own
`sweepCases/` subfolder. A tutorial-record sweep does not strip
`archive_dir_name` from its study before resolving study names (unlike a
factory-tutorial sweep, which does), so naming it here is refused as an
unknown key -- proven with a real `describe` preview, not assumed. This
study's own `driverFoam sweep-run` output is not archived per case today;
each case's raw output lives under the sweep's own `output_dir`.

**Also corrected 2026-09-26:** `PIMPLE.nNonOrthogonalCorrectors` is absent
from this case's own `system/fvSolution` (only `nOuterCorrectors` is
committed). A tutorial record's direct study-key channel writes with
`add_if_missing=False` (proven against a real copy of this case's
`system/fvSolution`: attempting to write a missing key raises `KeyError`),
so the four `nonorth1`/`combined` variants of this study (and every
`correctorN80` case) preview cleanly but cannot be committed by a real
`sweep-run` until a tutorial record gains a way to add a key its native
document does not yet have. This is a real, open gap, not a Python-vocabulary
key to rewrite away -- see the omniD side's own report for the evidence.

## Tracking & Outputs

All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
