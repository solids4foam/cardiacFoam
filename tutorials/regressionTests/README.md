# regressionTests tutorial assets

This folder stores shared regression payloads for tutorial cases.

- `NiedererEtAl2012/` contains the slab regression overrides and checkpoint file.
- `singleCell/` contains the single-cell regression overrides and checkpoint file.

The active regression tests live with their tutorial cases as `regressionTest.sh`
scripts, matching the solids4foam tutorial regression pattern. This directory is
kept only for legacy regression fixtures and compatibility while tests are moved
into the tutorials themselves.
