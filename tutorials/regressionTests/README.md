# regressionTests tutorial assets

This folder stores shared regression payloads for tutorial cases.

- `NiedererEtAl2012/` contains the slab regression overrides and checkpoint file.
- `singleCell/` contains the single-cell regression overrides and checkpoint file.

Each runnable tutorial keeps a thin `runRegressionTest.sh` wrapper in its own case
folder, but the reference data now lives here so regression assets are separated
from the tutorial definitions themselves.
