# Idealized heart/torso ECG application controls

## Reused case and geometry

The maintained local and remote source is `feat/bath-transmission-mms` at
`f5d8011e`; its local checkout is `/Users/simaocastro/noFrontend_bathMMS`.
`tutorials/idealizedHeart/electroHeartBath` supplies the coupled case and
`../mesh/heartTorso` its hydrated shared inputs. The live older remote
`omnid/heart-in-bath` remains at `6b1ecb2f`. A cached historical
`github/feat/heart-in-bath` tracking ref is not the selected live branch.

The ECG branch reuses only the shared input bundle and adds
[the opt-in ECG example](../tutorials/idealizedHeart/electroHeartECG/README.md).
It does not merge the bath branch's solver changes or reactivate its coupled
source by default. Future bath-bidomain work can reuse this geometry.

The 299482 tetrahedra partition into myocardium 123617, blood 31192,
lungs 9680, chest wall 17027 and remaining torso 117966. The lungs are two
idealized cylindrical conductivity zones, not separately conformal surfaces.
The heart/torso interface is conformal. Native v2412 `checkMesh` passes;
four faces exceed 70 degrees nonorthogonality and remain disclosed in its log.
Zone partition and positive physical tensors are checked independently.
The physical tensor pair reproduces the supplied effective monodomain
conductivity to a maximum component difference of 1.8e-12 S/m.

Existing chest locations are preserved. RA/LA/LL are the documented torso
equivalents for absent limbs. Native sampling uses nearest cells for E2/E5
on the skin, at distances 2.33/3.37 mm. Both recording paths use the same
owning-cell observations and balanced Wilson weights.

## Implementation additions

This extension adds a utility and tutorial, with no new domain, propagation
solver or ECG-model class. It reuses `myocardiumDomainInterface::New`,
`electroVolumeFieldDomain`, `extracellularPotentialDomain`, `ecgDomain`,
`torsoECG` and `leadFieldData`; their existing responsibilities remain intact.

- `replayLeadFieldECG` creates an independent recovery-owned bidomain provider,
  loads stored support-mesh Vm through the production cell map, iterates the
  existing extracellular potential domain and invokes production torso sampling.
  It uses the same Gi/Ge, body field, heart region and Wilson weights as the
  reciprocal observer. No cardiac source step occurs. Internal/boundary source
  arrays and stored files are checked for preservation.
- The reciprocal builder's default face-neighbour method remains unchanged.
  The supplied heart has 122 rank-deficient myocardial face stencils. Explicit
  `myocardialGradientStencil pointFallback` uses distributed point-connected
  myocardial neighbours only in those cells, excluding bath values. The final
  required-rank guard remains active and provenance identifies the choice.
- `-checkMyocardialGradient` exercises an affine field on the actual mesh.
  Serial and six-rank checks both give maximum error 2.88e-11. This checks
  reconstruction completeness/affine exactness, not nonlinear-gradient order.
- The forward/reciprocal control uses `distanceWeightedHarmonic`. The original
  bath case's `conormalHarmonic` configuration is separate. The new utility
  rejects an incompatible interface model, exterior forcing/patches, unmatched
  heart/body inputs, excessive linear residual, uncompleted recovery and source
  mutation. No exact discrete transpose of the corrected assembly is claimed.

## Mathematical comparison

For a frozen myocardial voltage, let the passive conductivity be
`K = Gi + Ge` in the heart and the prescribed body conductivity outside it.
The insulated passive forward problem has weak form

```text
integral_Omega grad(eta) . K grad(phi) dV
    = -integral_heart grad(eta) . Gi grad(Vm) dV.
```

Potential/current transmission is part of that conductor problem. A reference
fixes its additive constant. Let `D_j(phi)` denote the electrode observation
with +1 at Ej and -1/3 at RA/LA/LL. The balanced weights make `D_j` independent
of the potential gauge. A unit-current reciprocal solution satisfies

```text
integral_Omega grad(eta) . K grad(z_j) dV = D_j(eta).
```

For the symmetric physical tensors, reciprocity gives

```text
D_j(phi) = -integral_heart grad(z_j) . Gi grad(Vm) dV
         = integral_heart L_j . grad(Vm) dV,
L_j      = -transpose(Gi) grad(z_j).
```

Here `z_j` is normalized by reciprocal current. `L_j` has units m^-2,
so cell volume times `L_j . grad(Vm)` has units volts. Native electrodes are
owning-cell observations with matching cell current injection, rather than
a resolved finite-area skin contact model. The builder constructs z once;
the source run or replay only applies the volume projection thereafter.

This continuum identity motivates the comparison. Exact finite-grid equality
would additionally require reciprocal assembly/projection to be the adjoint
of the actual forward source operator. The corrected finite-volume solve and
myocardial cell-gradient quadrature here are independently assembled; that
additional relationship is not established. The measured difference therefore
needs operator/mesh investigation, rather than being attributed automatically
to propagation, physiology or an unconverged passive solve. Pseudo-ECG uses
its separate free-space kernel and is not a reference for this bounded problem.

The bulk expression for K does not prescribe interpolation of Gi+Ge across
the heart/torso interface. In the implemented forward solve, a global
extracellular/body assembly is combined with the mapped, sealed intracellular
heart assembly; its Vm source uses the corresponding heart operator. In the
builder, internal myocardial faces carry Gi+Ge, while heart/torso faces carry
the distance-weighted harmonic extracellular/body coefficient. The interface,
gradient corrections and projection must therefore be included explicitly in
the remaining discrete-operator comparison.

## Executed evidence

| Source dataset | Window | Samples | Six-lead relative RMS forward/reciprocal difference |
|---|---|---|---|
| Fresh monodomain, six ranks | 5–40 ms | 8 | 8.686% |
| Fresh heart-only bidomain, six ranks | 0.1–0.2 ms | 2 | 7.874% |
| Stored earlier coupled bath-bidomain source | 5–40 ms | 8 | 7.741% |

The difference is normalised by the corresponding forward cohort. It is not
an error against an exact solution. The fresh short bidomain run is a software
path control, not a QRS propagation study. The older coupled source has
byte-identical mesh geometry/connectivity/zones; its producer revision is not
certified here and no new coupled-propagation verification is inferred.
A fresh full-40-ms heart-only bidomain run was intentionally interrupted for
local execution cost before its first snapshot; that is not a completed result
or a numerical failure. The complete 40-ms recording window is exercised on
the fresh monodomain and separately on the archived coupled source.

Online/replay lead-field differences on fresh sources are at output rounding:
4.93e-12 V for monodomain and 4.80e-16 V for the short bidomain case. Independent
Python checks reconstruct Wilson references from raw torso electrode samples,
verify complete stored-time coverage and finite results, compare online/replay
projection, preserve source hashes and form Wilson-referenced pseudo traces.
Pseudo-ECG is a separate recording model and is not the bounded-torso truth.

Tightening the fixed-point tolerances from absolute 1e-11 V/relative 1e-8 to
1e-12 V/1e-9 changes forward samples by at most 3.63e-13 V in the monodomain
cohort and 9.39e-14 V in the archived source. The linear solver tolerance stays
1e-12. This selected fixed-point control does not explain the measured
8.686%/7.741% differences; it is not a bound on every algebraic error.

## Disposition and remaining numerical work

The maintained setup and frozen-source recording workflow execute, with an
identified finite-mesh discrepancy. Do not label the native forward and
reciprocal recordings equivalent, or silently impose an arbitrary pass bound.
Review the corrected conductor/source-gradient assemblies and assess mesh
sensitivity before making a stronger agreement claim. The nonorthogonal,
anisotropic, jump-containing case extends the prior scalar orthogonal column;
its result does not invalidate or inherit the manufactured projection orders.

Both changed/new utilities compile on native v2412 and Linux ARM64 v2312.
Case execution and serial/six-rank affine tests are native v2412. Four existing
serial/MPI lead-field configuration/recovery controls pass with zero skips.
The final matched-input-name guard is compiled on both platforms and checked
with a fresh positive replay and deliberate mismatched-body rejection. Full
PR CI, review and publication archive/tag release remain separate.

After the native application records, PR #59's `04e6bd84` compatibility update
was merged. It explicitly reads sub-mesh fvSchemes/fvSolution on versions
before 2406. The affected Linux v2312 library rebuild and independent
two-snapshot frozen-source replay pass, including the source/Wilson/online
checks. Maximum forward difference from the native two-snapshot replay is
1.62e-14 V. This is a separate integration control, not a Linux source
propagation study or a repeat of every earlier feature regression. An initial
build invocation omitted the backend resolver; its environment failure is
retained separately from the successful rebuild.

The local runtime evidence is retained under
`/tmp/cardiacFoam-lead-field-current-numerics/idealized-torso`. A separately
hashed local archive is retained at
`ElectrophysiologySolver/data/lead_field_applications/2026-10-14-idealized-torso`
in the manuscript repository. It is an application archive, not a new paper
convergence study or an immutable publication release. Full/tight waveforms
precede the final input-name guard/header normalization, with unchanged
numerical operators; the earlier binary hash and final-wrapper positive/negative
checks are recorded separately. Source hashes for those first datasets were
captured after their initial replay and before the tight repeat. The independent
final-wrapper positive fixture captures hashes before replay. No retrospective
file-preservation claim is made for an earlier un-hashed interval.
Source windows, initial/tight outputs, partial-run interruption and constructor
setup failures are retained rather than presented as successful results.
