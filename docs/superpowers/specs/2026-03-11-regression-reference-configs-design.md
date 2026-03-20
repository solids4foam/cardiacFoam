# Regression Reference Configs — Design Spec

**Date:** 2026-03-11

## Problem

Regression tests (`runRegressionTest.sh`) in `singleCell` and `NiedererEtAl2012` compare solver output against hardcoded reference values. Those values are only valid for a specific set of simulation parameters. If a user edits `constant/electroProperties` (e.g., changes `ionicModel`, `tissue`, or stimulus/physics params), the test runs with the modified config and produces different output, causing spurious failures.

## Goal

Ensure regression tests always run with the exact configuration the reference values were captured with, regardless of user edits.

## Design

### Reference config files

Each tutorial stores a locked copy of `constant/electroProperties` as `constant/electroProperties.reference`. This file is the authoritative config for the regression test. It is never edited by users — only updated when reference values are recaptured.

| Tutorial | Locked file |
|---|---|
| `tutorials/singleCell` | `constant/electroProperties.reference` |
| `tutorials/NiedererEtAl2012` | `constant/electroProperties.reference` |

### Restore function in `runRegressionTest.sh`

Each `runRegressionTest.sh` gains a `restoreReferenceConfigs` function that copies all `constant/*.reference` files to their canonical names before running the solver:

```bash
restoreReferenceConfigs() {
    for refFile in constant/*.reference; do
        [ -f "$refFile" ] || continue
        target="${refFile%.reference}"
        echo "Restoring ${target} from ${refFile}"
        cp "$refFile" "$target"
    done
}
```

This is called immediately before the solver invocation. The glob pattern means any future `.reference` file added to `constant/` is automatically picked up without further script changes.

### No solver changes required

The restore happens at the shell level before `runApplication`. The solver sees canonical configs as normal.

## Updating references

When reference values need to be recaptured (e.g., after a model change):

1. Update `constant/electroProperties` to the new canonical config
2. Run the solver and capture new output values
3. Update `system/singleCell.reference` or `system/smokeCheck.reference` with new values
4. Copy `constant/electroProperties` → `constant/electroProperties.reference`

## Files changed

- `tutorials/singleCell/constant/electroProperties.reference` — new file (copy of current `electroProperties`)
- `tutorials/singleCell/runRegressionTest.sh` — add `restoreReferenceConfigs` + call
- `tutorials/NiedererEtAl2012/constant/electroProperties.reference` — new file (copy of current `electroProperties`)
- `tutorials/NiedererEtAl2012/runRegressionTest.sh` — add `restoreReferenceConfigs` + call
