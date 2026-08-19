# Technical Debt: Orchestration Logic Duplication in Tutorial Specs

## Finding
During an architectural audit of the `driverFOAM` tutorial specs, we identified significant repetition of orchestration logic that should be centralized. 

Specifically, the function `_replace_blockmesh_resolution` is manually redefined and copy-pasted across 5 separate tutorial specifications:
1. `manufactured_eikonal_ecg.py`
2. `manufactured_fda.py`
3. `manufactured_fda_bath_bidomain.py`
4. `monodomain_and_eikonal_1d_cable_cv_convergence.py`
5. `niederer_2012.py`

In cases like the 1D Cable convergence and Niederer tutorials, the spec ignores the centralized `openfoam_driver.specs.utils.replace_single_block_mesh_resolution` entirely, opting to manually rewrite raw text-parsing loops (e.g., searching for `"hex (0 1 2 3 4 5 6 7)"`). 

## Implication
This violates `driverFOAM`'s centralized generation philosophy (where mesh provisioning math lives exclusively in `mesh_provisioning.py`). If the `blockMeshDict` format ever changes, or if `driverFOAM` introduces a new spatial tutorial, this duplicated logic will cause silent failures and maintenance overhead.

## Recommendation
Future agents or maintainers should:
1. Audit all tutorial specs for the `_replace_blockmesh_resolution` function.
2. Refactor them to strictly call `openfoam_driver.specs.utils.replace_single_block_mesh_resolution`.
3. Ensure any new spatial math (such as S1-S2 arrays) is written as a centralized utility in `openfoam_driver/specs/` before being adopted by a tutorial plugin.
