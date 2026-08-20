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

## Convergence Verification

Each case emits `setup/results/<key>_convergence.csv` (fresh, gitignored),
diffed within tolerance against the committed `reference/<key>_convergence.csv`,
plus a per-case `provenance.json` (OpenFOAM version + git SHAs + dict hashes).

| key | case |
|---|---|
| tet | monodomainPseudoECG (`setup/studies/tetConvergence`) |
| eikonal_tet | eikonalECG (`setup/studies/tetConvergence`) |
| coupling | monodomain1D3D |
| eikonal | eikonalECG |
| mono_spatial | monodomainPseudoECG |
| pseudo_ecg_spatial | monodomainPseudoECG |
| bidomain | bidomain |
| bidomain_tet | bidomain (`setup/studies/tetConvergence`) |
| bath | bathBidomain |
| bath_tet | bathBidomain (`setup/studies/tetConvergence`) |
| niederer | NiedererEtAl2011verification |

Cases whose sweep has not been run yet (`bidomain`, `bath`, and `bath_tet`
pending an interface-column confirmation) report **SKIP** until their
`reference/` CSV is frozen.
