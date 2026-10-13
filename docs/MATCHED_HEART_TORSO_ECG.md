# Matched heart/torso ECG comparison

## Purpose

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
a torso. The proposal uses upper-right/upper-left anterior torso sites for
RA/LA and a lower-left anterior site for LL. Each is 20 mm away from the
cranial/caudal and lateral truncation edges and 1 mm beneath the anterior
surface, permitting unambiguous interior point ownership.

| Reference | x (m) | y (m) | z (m) |
|---|---|---|---|
| RA | -0.119414800 | 0.149446125 | -0.104787473 |
| LA | -0.039891822 | -0.054769778 | -0.124073761 |
| LL | 0.183442648 | 0.017795990 | 0.028425169 |

These are declared modified torso placements, not certified full-limb
clinical electrode coordinates. VTK cell-location checks place all three
inside the existing 175,865-tetrahedron torso mesh. This verifies geometry;
production OpenFOAM/MPI point ownership and matched recordings remain to be
exercised when the tutorial variant is assembled.

The inspection, front/oblique renders, coordinate record and a balanced
`electrodePositions`/`leads`/`leadFields` fragment are retained locally under
`/tmp/cardiacFoam-lead-field-current-numerics/torso-placement`.
The fragment aliases the retained V1–V6 sites as E1–E6, producing output
channels V1–V6 with Wilson weights. It must be combined with the appropriate
physical conductor inputs; it does not define or run that conductor itself.

## Relation to the Niederer dataset

The supplied KCL/Niederer patient meshes have their own coordinate frames.
Their metadata arrays are organ tags and fibres; the accompanying information
PDF describes organ conductivities and ICD coils/cans, without a standard
chest/limb coordinate table. The idealized-case positions above belong to the
case's rotated box and must not be transplanted into a patient mesh without
registration and unit checks. An anatomical variant can preserve its own
trusted chest sites and make the same declared torso-reference assumption.
No private patient mesh or generated anatomical result is bundled here.
