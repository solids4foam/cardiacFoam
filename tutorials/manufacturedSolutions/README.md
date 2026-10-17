# Manufactured Solutions Tutorials

This folder groups the manufactured-solution verification cases by model scope.

## Cases

- `monodomainPseudoECG` : spatial manufactured-solution verification for the monodomain solver with pseudo-ECG output
- [monodomainLeadFieldECG](monodomainLeadFieldECG/README.md) : monodomain weighted-vector ECG functional with an independent analytical reference
- [bidomainLeadFieldECG](bidomainLeadFieldECG/README.md) : bidomain weighted-vector ECG functional with the same analytical reference
- `bidomain` : spatial manufactured-solution verification for the bidomain solver
- `bathBidomain` : spatial manufactured-solution verification for the bidomain solver with bath
- `eikonalECG` : activation-time and ECG manufactured verification for the eikonal solver
- `monodomain1D3D` : manufactured 1D-3D monodomain coupling verification
- `monodomainTotalLagrangianEM` : full coupled electromechanics manufactured-solution verification
- [`monodomainPseudoECG/sphericalShell`](monodomainPseudoECG/sphericalShell/README.md) : 3D monodomain with sealed curved walls crossed by a rotated full conductivity tensor, on a cubed sphere or gmsh tetrahedra, with the sealed-wall gradient and a parallel hex regression
- [`bidomain/sphericalShell`](bidomain/sphericalShell/README.md) : 3D heart-only bidomain with sealed curved walls crossed by two unequal rotated full conductivity tensors, on a cubed sphere or gmsh tetrahedra, with the sealed-wall gradient, a compatibility-projected phiE equation and a parallel hex regression
- [`eikonalECG/sphericalShell`](eikonalECG/sphericalShell/README.md) : 3D eikonal with sealed curved walls crossed by a rotated full conductivity tensor, on a cubed sphere or gmsh tetrahedra, with a pinned stimulus ball, the sealed-wall gradient and a parallel hex regression

## Naming Pattern

The subdirectory names follow the spatial PDE target being verified.

## Driver Usage

Driver commands below assume they are run from the repository root. The
repository-local wrapper is the supported entry point.

### Validate or run one case

Start with one normal strict entry to check the dictionaries and workflow for
the selected manufactured tutorial:

```bash
[omnidriver command to run]
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
described in omnidriver's own README, then check a sweep
without launching OpenFOAM. Spec: `tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/cartesianConvergence/sweep_hex_convergence.json`.

```bash
[omnidriver command to run]
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

Each case emits `setup/results/<key>_convergence.csv` (fresh, gitignored) plus a
per-case `provenance.json` (OpenFOAM version + git SHAs + dict hashes). There is
nothing committed to diff those against: the frozen `reference/<key>_convergence.csv`
baselines have been removed, and the harness that compared them never lived in this
repository. Automated coverage of these cases is the per-case `regression/`
directory, run through `tutorials/Alltest-regression`.

| key | case |
|---|---|
| tet | monodomainPseudoECG (`setup/studies/tetConvergence`) |
| eikonal_tet | eikonalECG (`-mesh tet`) |
| coupling | monodomain1D3D |
| eikonal | eikonalECG |
| mono_spatial | monodomainPseudoECG |
| pseudo_ecg_spatial | monodomainPseudoECG |
| bidomain | bidomain |
| bidomain_tet | bidomain (`-mesh tet`) |
| bath | bathBidomain |
| bath_tet | bathBidomain (`-mesh tet`) |
| niederer | NiedererEtAl2011verification |

## Bulk, boundary and gradient rows

Every field verifier (monodomain, bidomain, bath bidomain, eikonal) adds these
rows to its error summary for each potential `q` it reports, after its own rows:

| row | quantity |
|---|---|
| `q_bulk`, `q_boundary` | the error of `q` over the cells without, and with, a face on a physical patch of the mesh `q` lives on |
| `gradQ`, `gradQ_bulk`, `gradQ_boundary` | the error of the cell gradient of the solved `q`, with the solver's gradient scheme and boundary values |
| `gradReconQ`, `gradReconQ_bulk`, `gradReconQ_boundary` | the error of the cell gradient of the exact `q`, with exact cell and boundary-face values, through the same scheme: the reconstruction alone |

The physical patches are the sealed wall on a heart sub-mesh, the walls of a
whole-mesh case and the outer patches of a heart-and-bath mesh; processor,
cyclic and empty patches are not. L1 and L2 of every part are normalised by the
total volume of the mesh, so L1 and L2² of the parts add up to the total;
Linf is the maximum over the part. Gauge-shifted potentials (`phiE_gauge`)
keep their field rows under that name and their gradient rows under `phiE`.

On a heart-and-bath mesh the extracellular potential has a third part,
`phiE_interface`, the cells next to a heart–bath face. Its solved gradient is
the one the bath operator consumes (`phiEInterfaceGradient`: jump-aware under
`conormalHarmonic`), and the reconstruction rows leave the interface band out,
since the exact potential has a normal-derivative jump there.

Comparing `gradQ_boundary` with `gradReconQ_boundary` separates an
inconsistent boundary value of the solved field from the accuracy of the
reconstruction next to the boundary.
