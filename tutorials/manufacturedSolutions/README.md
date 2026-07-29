# Manufactured Solutions Tutorials

This folder groups the manufactured-solution verification cases by model scope.

## Cases

- `monodomainPseudoECG` : spatial manufactured-solution verification for the monodomain solver with pseudo-ECG output
- `bidomain` : spatial manufactured-solution verification for the bidomain solver
- `bathBidomain` : spatial manufactured-solution verification for the bidomain solver with bath
- `eikonalECG` : activation-time and ECG manufactured verification for the eikonal solver

## Naming Pattern

The subdirectory names follow the spatial PDE target being verified.

## Driver Usage

Run the whole manufactured suite:

```bash
./Allrun
```

Run only selected cases:

```bash
./Allrun monodomainPseudoECG bidomain
./Allrun bathBidomain
```

## Paper I reproducibility

Regenerate and check every reported convergence table with one command
(OpenFOAM sourced; needs only bash + git + stdlib `python3`):

```bash
./reproduce_paperI.sh                    # all cases (from repo root)
./reproduce_paperI.sh --skip-run tet     # re-check one case against its reference
```

Each case emits `studies/results/<key>_convergence.csv` (fresh, gitignored),
diffed within tolerance against the committed `reference/<key>_convergence.csv`,
plus a per-case `provenance.json` (OpenFOAM version + git SHAs + dict hashes).
The driver is data-driven by `applications/scripts/paperI_results/paperI_cases.tsv`.

| key | case | paper artefact |
|---|---|---|
| tet | monodomainPseudoECG (`studies/mesh/tet`) | tbl-tet-monodomain, tbl-tet-pseudoecg |
| eikonal_tet | eikonalECG (`studies/mesh/tet`) | §unstructured eikonal-tet |
| coupling | monodomain1D3D | coupled 1D–3D convergence |
| eikonal | eikonalECG | tbl-eikonal-*, tbl-eikonal-ecg-integral |
| mono_spatial | monodomainPseudoECG | tbl-monodomain-vm/-aux |
| pseudo_ecg_spatial | monodomainPseudoECG | tbl-pseudo-ecg |
| bidomain | bidomain | tbl-bidomain |
| bidomain_tet | bidomain (`studies/mesh/tet`) | tbl-tet-bidomain |
| bath | bathBidomain | tbl-bath-bidomain |
| bath_tet | bathBidomain (`studies/mesh/tet`) | tbl-bath-bidomain-tet |
| niederer | NiedererEtAl2011verification | fig-slab, tbl-niederer |

Cases whose sweep has not been run yet (`bidomain`, `bath`, and `bath_tet`
pending an interface-column confirmation) report **SKIP** until their
`reference/` CSV is frozen.
