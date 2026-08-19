# coupling - bathBidomain

## Purpose

Compares decoupled/baseline vs. predictor-corrector bath coupling
(`bidomainSolverCoeffs.bathPredictorCorrector`) at `N=10,20,40` on the
conformal tetrahedral mesh. Source of the paperI `bath_bidomain_tet_conformal`
experiment (`@tbl-bath-bidomain-corrector`-style sensitivity, not a spatial
convergence study).

Formerly `setup/mesh/tet/studies/coupling/run_coupling_study.sh` (relocated
here alongside its own study, matching this tutorial's other studies); mesh
generation and the tet electroProperties/fvSchemes overlay swap are handled
by driverFOAM's own `manufacturedBathBidomain` tet workflow DAG rather than
by hand-rolled bash.

## Execution

```bash
applications/scripts/driverFoam/bin/driverFoam sweep-run --spec setup/studies/coupling/sweep_coupling_study.json --output-dir <scratch_dir>
python3 setup/studies/coupling/summarize_coupling_study.py .
```

Run the summarizer from the tutorial's case root (`tutorials/manufacturedSolutions/bathBidomain`), not from `setup/studies/coupling/`. `--output-dir` for `sweep-run` is required by the CLI but only holds run-tracking state (`sweep_manifest.json`, per-case `run_document.json`) — the actual OpenFOAM data lands in the case root itself, described next.

### Where the output actually lands (verified, not the naive reading of `archive_dir_name`)

Each case runs **in-place** in the shared case root (`tutorials/manufacturedSolutions/bathBidomain/`), and the sweep engine archives it to `<case_root>/<caseId>/<archive_dir_name>/`, where `<caseId>` comes from the spec's `case_id_template` (`"<number_cells>_<bath_predictor_corrector>"`, e.g. `10_False`, `10_True`) and `<archive_dir_name>` is this spec's own `setup/studies/coupling/results/sweepCases`. So a real N=10 run leaves:

```
tutorials/manufacturedSolutions/bathBidomain/10_False/setup/studies/coupling/results/sweepCases/bathBidomainInterfaceMetrics.csv
tutorials/manufacturedSolutions/bathBidomain/10_True/setup/studies/coupling/results/sweepCases/bathBidomainInterfaceMetrics.csv
```

This is *not* the same as `applications/scripts/paperI_results/aggregate.py`'s `_sweep_cases_and_manifest()` helper, which assumes a single shared `setup/studies/<study>/results/sweepCases/` directory populated by a postprocess-consolidation step — that step is a driverFOAM stub as of this writing (`sweep-run`'s own output prints `"postprocess": {"status": "stub", ...}`), so nothing currently populates the shared location. `summarize_coupling_study.py` and `aggregate.py::_bath_tet()` were both rewritten to read the real per-`<caseId>` layout above instead.

`summarize_coupling_study.py` doesn't clean up `<caseId>` directories between runs — remove stale `N_False`/`N_True` dirs at the case root yourself before a fresh sweep if you don't want old data mixed into the summary.

## Status

Verified with a real `driverFoam sweep-run` at `N=10` (both `baseline` and `predictor`, 2026-08-19): both cases complete, and `summarize_coupling_study.py` correctly reads the archived output and reports a genuine physical difference (predictor-corrector coupling reduces `heartPhiE_L2`/`bathPhiE_L2` by roughly 50% at N=10 relative to baseline — sane and paper-consistent in direction). `N=20,40` unrun in this session but use the identical mechanism. `nOuterCorrectors 1`/`nNonOrthogonalCorrectors 1` are left at this case's own tet-overlay defaults rather than force-set (the old script explicitly forced both to 1; the checked-in default already matches, per the top-level README).

This also required a fix: the checked-in `constant/electroProperties` and `setup/mesh/tet/electroProperties` were both missing the `bidomainSolverCoeffs.{verificationModel,manufacturedBidomain}.fdaBathVariant` key that `_apply_case` always writes — every driverFOAM sweep for this tutorial (tet or hex, old specs included) crashed with a `KeyError` before this was added.

## Tracking & Outputs

All generated outputs are saved to the local `results/` folder, gitignored.
Do not commit generated OpenFOAM data.
