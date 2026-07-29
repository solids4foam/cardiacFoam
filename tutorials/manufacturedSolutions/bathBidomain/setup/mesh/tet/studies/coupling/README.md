# Bath PDE-coupling study

This study isolates one question in the conformal tetrahedral bath-bidomain
MMS without changing the ionic advance, physical timestep, or mesh-correction
policy:

1. Does enabling the predictor--corrector materially change the potential or
   interface-flux errors?

The tested controls are:

| Variant | `bathPredictorCorrector` | PDE sequence per physical step |
|---|---|---|
| `baseline` | `false` | update `phiE`, then solve `Vm` |
| `predictor` | `true` (default) | `Vm` predictor, `phiE`, `Vm` corrector |

The reaction/ionic model is advanced exactly once in every variant. Each
implicit `Vm` corrector uses the same old-time field, so the additional solves
do not advance physical time repeatedly. Both `phiE` and `Vm` use
`PIMPLE/nNonOrthogonalCorrectors`.

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
bash tutorials/manufacturedSolutions/bathBidomain/setup/mesh/tet/setup/coupling/run_coupling_study.sh
```

Useful overrides are `RESOLUTIONS="10 20 40"`,
`VARIANTS="baseline predictor"`, `RESULTS_DIR=/absolute/path`, and
`KEEP_WORK=1`.

Interpret potential and flux metrics separately. A coupling variant can
reduce the cell-centred potential error while leaving the local constitutive
face-flux error almost unchanged; algebraic single-valuedness of an assembled
internal-face flux is a separate conservation property.
