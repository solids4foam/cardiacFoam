# Bath PDE-coupling and non-orthogonal-correction study

This study isolates two questions in the conformal tetrahedral bath-bidomain
MMS without changing the ionic advance or physical timestep:

1. Does replacing the historical one-pass `phiE -> Vm` split with a stronger
   PDE coupling materially change the potential or interface-flux errors?
2. Does iterating the explicit correction of the global corrected Laplacian
   materially change those errors?

The tested controls are:

| Variant | PDE sequence per physical step | Global `phiE` non-orthogonal assemblies |
|---|---|---:|
| `baseline` | one `phiE -> Vm` pass | 1 |
| `predictor` | `Vm` predictor, `phiE`, `Vm` corrector | 1 |
| `phi8` | one `phiE -> Vm` pass | 9 |

The reaction/ionic model is advanced exactly once in every variant. Each
implicit `Vm` corrector uses the same old-time field, so the additional solves
do not advance physical time repeatedly.

The runner creates temporary cases, regenerates each tetrahedral mesh, forces
the paper formulation (`matchedSubmesh` plus
`distanceWeightedHarmonic`), and leaves the source tutorial untouched.

Before running, build the lightweight model (no `solids4foam` dependency):

```bash
source /Volumes/OpenFOAM-v2412/etc/bashrc
export FORCE_LIGHTWEIGHT_PHYSICSMODEL=1
source ./etc/resolveSolids4Foam.sh
wclean libso src/electroModels
wmake libso src/electroModels
```

Then run the default two-mesh screen:

```bash
bash tutorials/manufacturedSolutions/bathBidomain/setup/mesh/tet/studies/coupling/run_coupling_study.sh
```

Useful overrides are `RESOLUTIONS="10 20 40"`,
`VARIANTS="baseline predictor"`, `RESULTS_DIR=/absolute/path`, and
`KEEP_WORK=1`. The `phi8` variant is intentionally not recommended for the
large meshes until the inexpensive screen shows that it matters.

Interpret potential and flux metrics separately. A coupling variant can
reduce the cell-centred potential error while leaving the local constitutive
face-flux error almost unchanged; algebraic single-valuedness of an assembled
internal-face flux is a separate conservation property.
