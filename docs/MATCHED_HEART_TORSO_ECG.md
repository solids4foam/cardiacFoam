# Matched heart/torso ECG comparison

## Purpose

The protocol below is implemented in the opt-in
[`electroHeartECG` tutorial](../tutorials/idealizedHeart/electroHeartECG/README.md),
using the shared heart/torso inputs from `feat/bath-transmission-mms` at
`f5d8011e`. The [application report](IDEALIZED_HEART_ECG_VERIFICATION.md)
records executed controls and the remaining forward/reciprocal discrepancy.

Hold one heart/torso geometry, physical Gi/Ge and exterior conductivity,
electrode definitions, boundary/interface conditions and myocardial voltage
source fixed. Compare two ways of observing that source:

1. Solve the passive forward conductor problem for each stored voltage sample,
   then sample the resulting extracellular/torso potential at the electrodes.
2. Construct reciprocal sensitivities once for those same leads, then apply
   the existing volume integral to the same voltage-gradient samples.

For a Wilson-referenced chest channel, both recordings must use
`phi(chest_j) - (phi(RA)+phi(LA)+phi(LL))/3`. The reciprocal construction uses
+1 at the chest site and -1/3 at each reference site. Wilson's central terminal
is an average potential, not a spatial point to place in the mesh. See the
[AHA/ACC/HRS statement](https://doi.org/10.1016/j.jacc.2007.01.024).

Using the same stored Vm makes this an observation/conductor comparison.
Comparing an independently propagated monodomain source with a coupled
bath-bidomain source would also include source-model/bath-loading differences.
Those questions should be assessed separately. Pseudo-ECG may be recorded as
an additional free-space observation, but its different conductor model does
not supply the reference for this matched reciprocal comparison.

Assess complete traces, common-time RMS/maximum differences, gauge/reference
consistency, source preservation, mesh and reciprocal-solver sensitivity.
The existing static column is one bounded control. Its agreement does not
imply an exact discrete adjoint or full time-series equivalence on a larger
rotated/nonorthogonal mesh.

## Existing six case positions

The user's existing positions were found in the heart-in-bath branch,
`tutorials/idealizedHeart/electroHeartBath/constant/electroProperties`, at
`origin/omnid/heart-in-bath` revision 6b1ecb2f. They also occur in the saved
local run's `electroProperties.withDefaultValues`. They are retained exactly:

| Chest site | x (m) | y (m) | z (m) |
|---|---|---|---|
| V1 | 0.010407 | 0.056280 | -0.050468 |
| V2 | 0.030086 | 0.016628 | -0.061925 |
| V3 | 0.052524 | -0.019359 | -0.050439 |
| V4 | 0.070768 | -0.037057 | -0.026220 |
| V5 | 0.073474 | -0.043547 | -0.003609 |
| V6 | 0.068323 | -0.043011 | 0.018866 |

The case uses a rotated idealized torso box: its e1 axis points caudal, e2
patient-right and e3 posterior. The box is 0.32 x 0.26 x 0.20 m in that
frame. Its geometry was constructed around these chest sites; some are
slightly beneath its surface. Their existing positions are not altered.

The historical case samples point potentials at sites named V1–V6. Adding
Wilson referencing requires the three reference measurements or reciprocal
fields carrying their balanced weights; the point labels alone do not perform
that subtraction.

## Torso equivalents for absent limbs

The owner accepts assumed arm/leg references because the model contains only
a torso. The implemented case uses upper-right/upper-left anterior torso sites for
RA/LA and a lower-left anterior site for LL. Each is 20 mm away from the
cranial/caudal and lateral truncation edges and 1 mm beneath the anterior
surface, permitting unambiguous interior point ownership.

| Reference | x (m) | y (m) | z (m) |
|---|---|---|---|
| RA | -0.119414800 | 0.149446125 | -0.104787473 |
| LA | -0.039891822 | -0.054769778 | -0.124073761 |
| LL | 0.183442648 | 0.017795990 | 0.028425169 |

These are declared modified torso placements, not certified full-limb
clinical electrode coordinates. Earlier VTK inspection placed all three
inside the torso. Production OpenFOAM sampling on the complete 299482-cell
heart/torso mesh now exercises these references and the retained chest sites.
E2/E5 use bounded nearest-cell ownership at 2.33/3.37 mm from their supplied
skin points; other sites use interior ownership. The case declares a 15 mm
maximum distance, and both recording paths use the same electrode cells.

The inspection, front/oblique renders, coordinate record and a balanced
`electrodePositions`/`leads`/`leadFields` fragment are retained locally under
`/tmp/cardiacFoam-lead-field-current-numerics/torso-placement`.
The original fragment aliases the retained V1–V6 sites as E1–E6, producing
output channels V1–V6 with Wilson weights. Those definitions and the physical
conductor inputs now ship together in the maintained tutorial.

## Executed matched recordings

Fresh monodomain sources supply eight saved samples from 5 to 40 ms, with
both online observers. A fresh heart-only bidomain software-path check supplies
two samples at 0.1/0.2 ms; an earlier stored coupled bath-bidomain source supplies
another eight-sample 5–40 ms recording control. Its producer revision is not
certified here. These are three distinct datasets, not a propagation comparison.

Independent checks reconstruct Wilson references from raw torso samples and
compare online/replayed reciprocal projection. The six-lead relative RMS
forward/reciprocal differences are respectively 8.686%, 7.874% and 7.741%.
Tighter recovery fixed-point controls barely change the forward traces and
do not explain these differences. The complete workflow executes, but a review
of the corrected conductor/source-gradient assemblies and mesh sensitivity
remains necessary before claiming stronger recording agreement.

## Relation to the Niederer dataset

The supplied KCL/Niederer patient meshes have their own coordinate frames.
Their metadata arrays are organ tags and fibres; the accompanying information
PDF describes organ conductivities and ICD coils/cans, without a standard
chest/limb coordinate table. The idealized-case positions above belong to the
case's rotated box and must not be transplanted into a patient mesh without
registration and unit checks. An anatomical variant can preserve its own
trusted chest sites and make the same declared torso-reference assumption.
No private patient mesh or generated anatomical result is bundled here.
