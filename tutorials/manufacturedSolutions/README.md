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
