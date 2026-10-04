# coupling - bathBidomain

## Purpose

Compares decoupled/baseline vs. predictor-corrector bath coupling
(`bidomainSolverCoeffs.bathPredictorCorrector`) at `N=10,20,40,80` on the
conformal tetrahedral mesh. Source of the paperI `bath_bidomain_tet_conformal`
experiment (`@tbl-bath-bidomain-corrector`-style sensitivity, not a spatial
convergence study).

Mesh generation is the omniD `manufacturedBathBidomain` record's tet route (see
the top-level README).

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

The study runs through the `manufacturedBathBidomain` tutorial record, which
stages every case under the sweep's `--output-dir` as `cases/<caseId>/` (e.g.
`cases/10_False/postProcessing/bathBidomainInterfaceMetrics.csv`). Point
`summarize_coupling_study.py` at the sweep's `cases/` directory.

## Status

Only `N=10` (both `baseline` and `predictor`) has been run; `N=20,40,80` use the
same mechanism. The `N=80` pair is needed to tell a bath-coupling sensitivity
from the separate finest-level interface-current anomaly. `nOuterCorrectors 1`/
`nNonOrthogonalCorrectors 1` are left at the case's own defaults (see the
top-level README).

The study sets `bidomainSolverCoeffs.bathPredictorCorrector` directly
(`false`/`true`; the case's own value is `yes`), and the `groundElectrode`
variant as the top-level README shows.

## Tracking & Outputs

All generated outputs are saved to the local `results/` folder, gitignored.
Do not commit generated OpenFOAM data.
