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

```bash
[omnidriver command to run]
python3 tutorials/manufacturedSolutions/bathBidomain/setup/studies/coupling/summarize_coupling_study.py <sweep output dir>
```

Run the summarizer from the repository root. It writes `raw_results.csv` and
`summary.md` into the sweep output directory.

### Where the output actually lands

The study runs through the `manufacturedBathBidomain` tutorial record, which
stages every case under the sweep's output directory as `cases/<caseId>/` (e.g.
`cases/10_False/postProcessing/bathBidomainInterfaceMetrics.csv`), with the
case's axis values in `<caseId>/case_record.json`. Point
`summarize_coupling_study.py` at the sweep's output directory.

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

Generated outputs stay in the sweep output directory.
Do not commit generated OpenFOAM data.
