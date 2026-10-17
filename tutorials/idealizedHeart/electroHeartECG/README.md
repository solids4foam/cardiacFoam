# Matched ECG recordings on the idealized heart and torso

One shared heart/torso mesh supports monodomain and heart-only bidomain
source runs, concurrent reciprocal and pseudo-ECG recording, and subsequent
passive torso recovery on exactly the saved Vm snapshots. The passive
conductor never feeds back into either source run.

## Shared geometry and provenance

The supplied `../mesh/heartTorso` bundle comes from
`feat/bath-transmission-mms` revision `f5d8011e`, whose maintained local checkout
is `noFrontend_bathMMS`. The older `omnid/heart-in-bath` branch is the initial
version; the selected bundle includes the later sealed-heart and conductivity
updates. Its heart cells and canonical heart/anatomy inputs are shared with
`electroHeart`. The future coupled bath-bidomain tutorial can use the same
geometry without requiring a different torso.

| Region | Cells | Passive conductivity (S/m) |
|---|---|---|
| Myocardium | 123617 | Physical intracellular/extracellular tensor fields |
| Blood | 31192 | 0.7 |
| Lungs | 9680 | 0.0389 |
| Chest wall | 17027 | 0.05 |
| Remaining torso | 117966 | 0.22 |

The five zones partition 299482 tetrahedra. Lungs are two idealized cylindrical
conductivity zones on that mesh; their surfaces are not separately conformal.
The heart/torso interface is conformal. The scalar body-field value in heart
cells is ignored: physical Gi/Ge supply the heart operator. The supplied
monodomain tensor equals the series combination of the commuting physical
Gi/Ge tensors to input precision.

The existing E1–E6 chest positions are retained. RA/LA/LL are declared anterior
torso equivalents for missing limbs, not patient-specific limb coordinates.
`constant/electrodes` records them. Each channel is
`Vj = phi(Ej) - (phi(RA)+phi(LA)+phi(LL))/3`; Wilson's reference is a virtual
average potential. The same balanced weights drive reciprocal construction
and combine the forward torso samples. The explicit 15 mm maximum distance
permits bounded nearest-cell placement for a chest point on the skin;
interior point ownership is used when available.

## Source and conductor variants

Both source variants retain the baseline Bueno–Orovio ionic heterogeneity,
Purkinje network and pre-pacing configuration. They solve only the myocardial
cellZone and use a sealed heart boundary. `monodomain` uses its effective
propagation tensor; `bidomain` uses the supplied physical pair and a reference
inside the heart. Their propagation is not assumed identical.

The reciprocal and replay conductor both use
`distanceWeightedHarmonic` interface interpolation. The existing coupled
bath tutorial selects `conormalHarmonic`; that is a separate configuration,
not silently substituted here. No exact finite-grid transpose relationship
between corrected forward recovery and cell-gradient projection is assumed.
Comparisons must report their measured discrepancy and solve/mesh sensitivity.

The canonical tetrahedral heart contains boundary cells without enough
myocardial face-neighbours for a 3-D reciprocal gradient. This case explicitly
selects `myocardialGradientStencil pointFallback`. Only rank-deficient cells
use additional point-connected myocardial neighbours; bath samples are
excluded, full-rank face stencils retain their existing reconstruction, and
the required-rank guard remains active. Native distributed stencils support
processor exchange; transformed coupled patches are rejected by that native
stencil implementation. The setup checks an affine field on the actual mesh
before building physical lead fields. This affine control is not a general
nonlinear-gradient accuracy or convergence-order claim.

## Run

Hydrate the shared Git-LFS assets and build the normal solver/utilities first.
Each fresh source run requires its own case directory, or `Allclean` before
switching variants.

```bash
./Allrun prepareOnly             # mesh check, affine check and six lead fields
./Allrun monodomain parallel     # 40 ms source, both online observers, replay
./Allrun replayOnly              # replay existing reconstructed Vm snapshots
./Allclean
./Allrun bidomain parallel -endTime 0.0002 -writeInterval 0.0001
# Omit the two time overrides for the default 40 ms bidomain source run.
```

`replayOnly` replaces the replay outputs and its log; preserve a comparison
dataset before repeating. It retains the original source-hash record and
rejects a changed source.

The configured source window is 0–40 ms with a 20 microsecond timestep and saved
snapshots every 5 ms. Explicit `-endTime` and `-writeInterval` options write
only `local/runControls`; a short execution check does not establish QRS
propagation accuracy or complete-window recording. Full-precision writes preserve mesh provenance through
decomposition and source values through replay. Parallel source runs are
reconstructed before the serial passive replay. `Allclean` removes generated
local inputs/results and leaves the canonical shared bundle untouched.

## Assess the recordings

`replayLeadFieldECG` reads `constant/passiveECGProperties`, creates an independent
recovery provider, and holds each saved myocardial Vm fixed. It records native
forward torso potentials combined with the Wilson weights and the reciprocal
projection of that same source. It rejects unconverged recovery, excessive
linear residual, incompatible exterior forcing/patches and source mutation.

`postProcessing/matchedLeadFieldECG.dat` contains every paired sample,
correction count, update, linear residual and source-preservation result.
`matchedLeadFieldECGSummary.dat` gives per-channel absolute RMS differences,
forward RMS amplitudes and maximum differences. Online `leadFieldECG` and
`pseudoECG` outputs remain separate; pseudo-ECG is a distinct free-space
observation and its raw electrode values require the same Wilson subtraction
before comparing lead waveforms.

Compare the online and replay reciprocal values at common stored times, then
assess the forward/reciprocal discrepancy and tighten recovery controls where
needed. A small discrepancy on this fixed mesh is supporting application
control, not an asymptotic-order, anatomical-validation or general discrete-
adjoint result. Independently propagated mono and coupled bath signals also
include source/bath-loading differences and do not isolate recording error.

The maintained checker captures source-file hashes before replay, checks them
again afterwards, independently forms Wilson leads from raw torso potentials,
checks source-preservation/linear-convergence records, and compares online and
replayed reciprocal values at common times. It reports absolute/per-channel
and cohort relative RMS differences rather than treating the forward operator
as an exact solution. It also writes `pseudoWilsonECG.dat`, applying the same
reference to the independent pseudo-ECG electrode samples.

After a source run, `./regression/regressionTest.sh --check-only` validates its
existing outputs through the shared regression framework. Without
`--check-only` it starts a fresh default monodomain run. This larger application
is opt-in rather than part of the default lightweight sweep. Its contract
checks do not impose an invented exact-accuracy bound on the measured
forward/reciprocal discrepancy.

The [executed application report](../../../docs/IDEALIZED_HEART_ECG_VERIFICATION.md)
records the fresh monodomain and short bidomain checks, the separately
qualified archived-source replay, and the open finite-mesh recording
discrepancy. Execution-contract passes do not establish ECG equivalence.
