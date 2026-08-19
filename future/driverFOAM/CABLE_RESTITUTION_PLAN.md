# Deterministic Implementation Plan: Spatial Restitution driverFOAM Spec

This document details the exact, deterministic steps to migrate the bespoke 1D cable S1-S2 protocol into a normalized `driverFOAM` tutorial specification.

## 1. Directory Cleanup
Execute the following commands from the repository root:
```bash
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/run_smart_restitution.sh
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/run_restitution.sh
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/run_calib.sh
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/run_phase2.sh
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/run_phase2_quick.sh
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/debug_cross.py
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/update_stimuli.py
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/print_vm.py
rm tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/log.*
```

## 2. Create Centralized Spatial Math
**File:** `applications/scripts/driverFoam/openfoam_driver/specs/spatial_pacing.py`

**Code to inject:**
```python
from __future__ import annotations

def generate_spatial_s1_s2_stimulus_lists(
    s1_interval_ms: float, n_s1: int, s2_interval_ms: float, n_s2: int,
    bounds_min: str, bounds_max: str, duration_s: str, intensity: str
) -> dict[str, str]:
    times = []
    for i in range(n_s1):
        times.append(i * (s1_interval_ms / 1000.0))
    for i in range(n_s2):
        times.append((n_s1 * (s1_interval_ms / 1000.0)) + i * (s2_interval_ms / 1000.0))
    
    count = len(times)
    return {
        "stimulusStartTimeList": "(" + " ".join(f"{t:.6g}" for t in times) + ")",
        "stimulusLocationMinList": "(" + " ".join([bounds_min] * count) + ")",
        "stimulusLocationMaxList": "(" + " ".join([bounds_max] * count) + ")",
        "stimulusDurationList": "(" + " ".join([duration_s] * count) + ")",
        "stimulusIntensityList": "(" + " ".join([intensity] * count) + ")",
    }
```

## 3. Create driverFOAM Tutorial Spec
**File:** `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/cable_restitution_curves.py`

**Code to inject:** (Standard spec boilerplate omitted for brevity; exact overrides shown)
```python
# Inside _apply_case:
from openfoam_driver.specs.utils import replace_single_block_mesh_resolution, set_end_time
from openfoam_driver.specs.spatial_pacing import generate_spatial_s1_s2_stimulus_lists
from openfoam_driver.plugins.cardiacfoam.overrides import apply_electro_property_overrides

def _apply_case(case_root, case, ...):
    # 1. Mesh
    replace_single_block_mesh_resolution(
        case_root / block_mesh_dict_relpath,
        float(case.params["dx_mm"]),
        "cable", # dimension
        resolution_by_dimension={"cable": defaults.CROSS_SECTION_CELL_COUNTS}
    )
    
    # 2. Pacing arrays
    stimulus_arrays = generate_spatial_s1_s2_stimulus_lists(
        s1_interval_ms=case.params["s1Interval"], n_s1=case.params["nS1"],
        s2_interval_ms=case.params["s2Interval"], n_s2=case.params["nS2"],
        bounds_min="(0 0 0)", bounds_max="(2e-3 2e-4 2e-4)",
        duration_s="4e-3", intensity="50000"
    )
    
    # 3. Apply overrides to electroProperties
    overrides = {
        f"monodomainSolverCoeffs.externalStimulus.{k}": v 
        for k, v in stimulus_arrays.items()
    }
    overrides["monodomainSolverCoeffs.ionicModel"] = str(case.params["ionicModel"])
    overrides["monodomainSolverCoeffs.tissue"] = str(case.params["tissue"])
    apply_electro_property_overrides(case_root / electro_properties_relpath, overrides)
    
    # 4. End time
    end_time = (case.params["s1Interval"] * case.params["nS1"] + case.params["s2Interval"] * case.params["nS2"]) / 1000.0 + 0.1
    set_end_time(case_root / control_dict_relpath, end_time)
```

## 4. Register Tutorial
**File:** `applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/registry.py`

**Modifications:**
- **Add import:** `from openfoam_driver.plugins.cardiacfoam.tutorials.cable_restitution_curves import make_spec as make_cable_restitution_curves_spec`
- **Add to `SPEC_FACTORIES` map:** `"cableRestitutionCurves": make_cable_restitution_curves_spec,`
- **Add to `REGISTERED_TUTORIALS` tuple:** `"cableRestitutionCurves",`

## 5. Normalized Post-Processing
**File:** `tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/setup/postProcessing_cableRestitution.py`

**Code to inject:**
```python
from pathlib import Path
import json

def run_postprocessing(output_dir: Path, **kwargs) -> None:
    """
    NOTE: driverFOAM's PostprocessTask/run_postprocess_tasks hook was removed
    (2026-08-18) after auditing found it unreachable from any live CLI action
    -- nothing in the run --strict / sweep-run execution engine ever called
    it. There is currently no automated post-DAG hand-off in driverFOAM; this
    function's signature is kept as a reasonable shape for whatever hand-off
    mechanism replaces it, but wiring it up will need a different integration
    point than "driverFOAM invokes it automatically."
    1. Parse postProcessing/cableProbes/*/Vm
    2. Identify zero-crossings
    3. Compute S2 CV
    4. Write outputs matching extract_cv.py format into `output_dir`.
    """
    # (Extract logic strictly ported from debug_cross.py)
    ...
```
