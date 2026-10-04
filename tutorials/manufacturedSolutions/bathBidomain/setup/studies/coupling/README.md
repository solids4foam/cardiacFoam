# coupling - bathBidomain

## Purpose

Compares decoupled/baseline vs. predictor-corrector bath coupling
(`bidomainSolverCoeffs.bathPredictorCorrector`) at `N=10,20,40,80` on the
conformal tetrahedral mesh. Source of the paperI `bath_bidomain_tet_conformal`
experiment (`@tbl-bath-bidomain-corrector`-style sensitivity, not a spatial
convergence study).

Formerly `setup/studies/tetConvergence/studies/coupling/run_coupling_study.sh` (relocated
here alongside its own study, matching this tutorial's other studies); mesh
generation is the omniD `manufacturedBathBidomain` record's tet route (see
the top-level README) rather than hand-rolled bash.

## Execution

Spec: `tutorials/manufacturedSolutions/bathBidomain/setup/studies/coupling/sweep_coupling_study.json`.

```bash
omnidriver --plugin cardiacfoam sweep-run --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/coupling/sweep_coupling_study.json --output-dir <scratch output dir>
python3 tutorials/manufacturedSolutions/bathBidomain/setup/studies/coupling/summarize_coupling_study.py tutorials/manufacturedSolutions/bathBidomain
```

`omnidriver` is the external orchestration add-on (not part of this repo;
see the root `CLAUDE.md`). Run the summarizer from the repository root.
`--output-dir` holds run-tracking state (`sweep_manifest.json`, per-case
`run_document.json`) while the actual OpenFOAM data lands in the case root
itself, described next.

### Where the output actually lands

Corrected 2026-09-26 (tutorials-are-pointers 5.4a): this section described
an in-place run archived under `<case_root>/<caseId>/<archive_dir_name>/`.
The study now runs through the `manufacturedBathBidomain` tutorial record,
which stages every case under the sweep's own `--output-dir`, as
`cases/<caseId>/` (e.g. `cases/10_False/postProcessing/bathBidomainInterfaceMetrics.csv`);
`archive_dir_name` is not a record study key and is dropped.
`summarize_coupling_study.py` still reads the old per-`<caseId>` layout under
the case root, so point it at the sweep's `cases/` directory.

## Status

`N=10` (both `baseline` and `predictor`) runs to completion; `summarize_coupling_study.py`
correctly reads the archived output and reports a genuine physical
difference (predictor-corrector coupling reduces `heartPhiE_L2`/`bathPhiE_L2`
by roughly 50% at N=10 relative to baseline — sane and paper-consistent in
direction). `N=20,40,80` use the identical mechanism and have not been run.
The `N=80` pair is required to distinguish a bath-coupling sensitivity from
the separate finest-level interface-current anomaly. `nOuterCorrectors 1`/
`nNonOrthogonalCorrectors 1` are left at this case's own tet-overlay
defaults rather than force-set, matching the checked-in default (see the
top-level README).

The study sets `bidomainSolverCoeffs.bathPredictorCorrector` directly
(`false`/`true`; the case's own value is `yes`), and the `groundElectrode`
variant as the top-level README shows.

## Tracking & Outputs

All generated outputs are saved to the local `results/` folder, gitignored.
Do not commit generated OpenFOAM data.
