# Artifact Prediction Bug in driverFOAM

## Issue Description
When `driverFOAM`'s strict planning component predicts expected artifacts for a tutorial, it calls the `cardiacFoam` plugin's `predict_cardiac_artifacts`.
This function attempts to detect the ionic and active tension export lists by reading the `electroProperties` template. 
It uses `detect_ionic_export_list` and `detect_active_tension_export_list` in `openfoam_driver/plugins/cardiacfoam/detection.py`.

The regex used extracts the variables inside the `export ( ... )` list. If the list is intentionally empty (e.g., `export ();`), the tokens array `tokens` becomes `()`. 
However, the functions return `tokens if tokens else None`. Since `()` is falsy, the functions return `None`.

When `detect_ionic_export_list` returns `None`, `predict_cardiac_artifacts` assumes that the export list could not be parsed or wasn't provided in the dictionary. It then falls back to the `ionic_model_catalog.py` and `active_tension_catalog.py` and appends all of the "recommended" exports for the given ionic model.

For example, if the model is `Stewart` and the `electroProperties` sets `export ();`, the predictor still injects variables like `Cai`, `AV_i_Na`, etc., because the function returned `None`. This subsequently causes the strict plan runner to fail with `missing_expected_artifacts` because OpenFOAM correctly produced zero output files for the `()` export.

## Solution
Modify `openfoam_driver/plugins/cardiacfoam/detection.py` to correctly return the empty tuple instead of `None` when the `export` block is explicitly empty:

```python
def detect_ionic_export_list(electro_properties_path: Path) -> tuple[str, ...] | None:
    # ...
    tokens = tuple(t for t in match.group(1).split() if t)
    return tokens  # <-- DO NOT return `tokens if tokens else None`
```

Do the exact same thing for `detect_active_tension_export_list`.

## Temporary Workaround
To bypass this bug without modifying `driverFOAM`'s core code, the tutorial's `electroProperties` template must provide at least one variable in the `export` list. For example, changing `export ();` to `export (Vm);` prevents the function from returning `None` and correctly overrides the catalog fallback.
