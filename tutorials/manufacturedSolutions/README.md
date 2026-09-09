# Manufactured Solutions Tutorials

This folder groups the manufactured-solution verification cases by model scope.

## Cases

- `monodomainPseudoECG` : spatial manufactured-solution verification for the monodomain solver with pseudo-ECG output
- `bidomain` : spatial manufactured-solution verification for the bidomain solver
- `bathBidomain` : spatial manufactured-solution verification for the bidomain solver with bath
- `eikonalECG` : activation-time and ECG manufactured verification for the eikonal solver
- `monodomain1D3D` : manufactured 1D-3D monodomain coupling verification
- `monodomainTotalLagrangianEM` : full coupled electromechanics manufactured-solution verification

## Naming Pattern

The subdirectory names follow the spatial PDE target being verified.

## Driver Usage

Driver commands below assume they are run from the repository root. The
repository-local wrapper is the supported entry point.

### Validate or run one case

Start with one normal strict entry to check the dictionaries and workflow for
the selected manufactured tutorial:

```bash
driverFoam plan \
    --strict \
    --entry manufacturedMonodomainPseudoECG
driverFoam run \
    --strict \
    --entry manufacturedMonodomainPseudoECG
```

`plan` validates and reports the generated case. `run` executes that one case
and uses the entry's configured case output location. Replace the entry with
`manufacturedBidomain`, `manufacturedBathBidomain`, `manufacturedEikonalECG`,
`manufacturedMonodomain1D3D`, or
`manufacturedMonodomainTotalLagrangianEM` as appropriate.

### Run a manufactured convergence sweep

Use a sweep only for a study that varies mesh, timestep, numerical scheme, or
another declared axis. Give every simultaneous run its own output directory;
the directory contains the manifest and per-case state, so reusing one while
another run is active can mix or resume state from the other run.

Before running a case, configure the host-specific OpenFOAM/runtime file as
described in the driverFOAM add-on's own README, then check a sweep
without launching OpenFOAM:

```bash
export DRIVERFOAM_RUNTIME_CONFIG=/absolute/path/driverfoam-runtime.yaml
driverFoam sweep-plan \
    --spec tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json \
    --output-dir .tmp/driverfoam/monodomainPseudoECG-cartesian
```

Use `sweep-run` only after the plan is ready. The current checkout was audited
on 2026-08-26: valid cases materialise, but strict execution is presently
blocked by local runtime preflight because a built library hash no longer
matches the `cardiacFoam` build manifest and the executable is older than the
current sources. Rebuild the OpenFOAM/plugin libraries and refresh the build
manifest before attempting a real solve. Several specs also have independent,
pre-OpenFOAM issues documented in their case READMEs; `sweep-plan` is the
required first diagnostic.

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
