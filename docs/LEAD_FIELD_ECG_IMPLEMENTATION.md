# Lead-field ECG: implementation and verification

This document describes the monodomain/bidomain feature in
`codex/lead-field-ecg-integration`, based on the sealed-wall implementation in
PR #59. It explains the runtime classes, preprocessing, numerical contract,
tutorials and the evidence required for the separate ECG pull request.

## What the feature computes

For a fixed conductor and a balanced recording channel, preprocessing solves
one reciprocal current problem and stores its myocardial sensitivity. During
cardiac simulation the observer evaluates

```text
L_l = -Gi^T grad(z_l)
y_l(t) = sum_cells volume * L_l dot grad(Vm(t))
```

The reciprocal potential `z_l` is normalized per ampere and has units ohm.
`Gi` and `Ge` are physical conductivities in S/m. The weighted vector `L_l`
has units m^-2 and excludes cell volume. The projection therefore produces
volts. Conductivity and current normalization are already in the saved vector:
applying Gi or a pseudo-ECG prefactor again would change the recording model.

For Wilson-referenced V1, the reciprocal electrode weights are +1 at E1 and
-1/3 at RA, LA and LL. V2-V6 use the same reference with E2-E6. The importer
reads named channel fields; it does not pretend these combinations are single
point electrodes. All six combinations are supplied in
[the Wilson example](../applications/utilities/buildLeadFields/examples/wilsonLeads).
Their coordinates must be provided for the actual conductor.

The passive continuum identity pairs the cardiac source with the reciprocal
conductor solution. It assumes compatible conductor/interface conditions and
no independent extracellular forcing. Preprocessing uses Gi+Ge in myocardium,
passive conductivity outside, and the existing extracellular-only face helper
at heart/passive interfaces. The numerical construction is a corrected
finite-volume solve plus myocardial reconstruction and cell quadrature. It is
not asserted to be the exact transpose of the full coupled bath matrix.

A torso mesh is needed to construct a torso-sensitive field. After explicit
myocardial export, the cardiac simulation can use only the heart and saved
vectors. This removes a surrounding recording solve from each timestep; it
does not recover bath-loading effects omitted from the source model.

## Why the observer belongs to ECG

The online operation reads Vm and does not advance another state or feed back
to ionic/PDE evolution. The existing `ecgSolver` selection and
`electroStateProvider` interface provide that contract for both source models.
A new potential-domain subclass is unnecessary for this read-only projection.

The existing potential-domain machinery is still used where an extracellular
potential really is solved, including the independent recovery test below.
No monodomain or bidomain PDE kernel was added or replaced by this feature.

## Added classes and responsibilities

| Class | Source | Responsibility |
|---|---|---|
| `leadFieldECG : ecgSolver` | [leadFieldECG.H/C](../src/electroModels/ecgModels/leadFieldECG/leadFieldECG.H) | Runtime selection, static initialization, channel names and projection delegation. |
| `leadFieldData` | [leadFieldData.H/C](../src/electroModels/ecgModels/leadFieldECG/leadFieldData.H) | Physical input loading, native vector import, addressing, provenance checks and volume projection. |
| `manufacturedLeadFieldECGVerifier : ecgVerificationModel` | [verifier H/C](../src/verificationModels/ecgVerification/manufacturedLeadFieldECGVerifier.H) | Independent analytical functional plus source/gradient diagnostics. |
| `affineState : electroStateProvider` | [fixture utility](../applications/utilities/setManufacturedLeadFields/setManufacturedLeadFields.C) | Verification-only adapter for a prescribed voltage, conductivity and optional extracellular snapshot. |

`leadFieldData` is a helper rather than another runtime solver. This lets the
production observer and verification utilities use the same import/projection
contract without duplicating its numerical implementation. `affineState` is
local to a test utility, not another cardiac model.

## Added applications

| Application | Location | Purpose |
|---|---|---|
| `buildLeadFields` | [utility](../applications/utilities/buildLeadFields/README.md) | Static passive-conductor solve, sensitivity generation, optional reciprocal potentials and provenance. |
| `setManufacturedLeadFields` | [utility](../applications/utilities/setManufacturedLeadFields/README.md) | Synthetic sensitivity fixture and private-gradient/source-preservation controls. |
| `checkLeadFieldRecovery` | [utility](../applications/utilities/checkLeadFieldRecovery/README.md) | Frozen-voltage production recovery and actual torso sampling against independent and reciprocal references. |

Each application has `Make/files`, `Make/options` and a utility manifest for the
existing discovery/build mechanism. They are command-line applications, not
new production-domain classes.

Within myocardium, the reciprocal gradient uses only myocardial neighbours,
including processor neighbours. It excludes bath-side slopes without imposing
an artificial zero current at the restriction. The normal matrix is ranked
with SVD using relative cutoff 1e-12. Full-rank stencils use direct inversion;
lower-dimensional stencils use the SVD range pseudoinverse. Unresolved active
geometric directions cause rejection. This replaces an unreliable cubic
eigenbasis reconstruction for rotated repeated eigenvalues.

## Existing files extended

- [ecgDomain.H/C](../src/electroModels/electroDomains/ecgDomain/ecgDomain.H):
  separates source electrodes from output channels, initializes the observer,
  routes optional extracellular recovery and manages sampling/restart output.
- [ecgSolver.H/C](../src/electroModels/electroDomains/ecgDomain/ecgSolver.H):
  provides initialization/channel hooks while preserving legacy defaults.
- [electrophysicsSystemBuilder.C](../src/electroModels/core/system/electrophysicsSystemBuilder.C):
  delegates observer/potential setup to the owning domain rather than
  duplicating that lifecycle in the factory builder.
- [conductivityFieldIO.H/C](../src/genericWriter/conductivityFieldIO.H):
  supports explicit physical Gi/Ge through the established field/uniform
  contract and its legacy-name handling.
- [ecgModelIO.H/C](../src/genericWriter/ecgModelIO.H):
  validates a restart trace before append: exact header, finite values and
  increasing times. Equality in the saved folder/trace representation permits
  rounding of the same step; genuinely later samples remain a rejection.
- [torsoECG.C](../src/electroModels/ecgModels/torsoECG/torsoECG.C):
  obtains cell centres collectively before owner-only sampling branches,
  avoiding lazy coupled geometry communication on only one MPI rank.
- Library `Make/files/options`, utility/model READMEs, the architecture record
  and tutorial/configuration indexes register and document the feature.

## Configuration and source preservation

```foam
ecgDomains
{
    ECG { ecgSolver leadFieldECG; }
    Pseudo { ecgSolver pseudoECG; sigmaExtracellular 0; }
}
```

The second observer is optional. Different solver types have different output
files. Use one observer per solver type for its default trace filename; add
more lead channels to the one `leadFieldProperties` configuration.

`constant/leadFieldProperties` specifies physical inputs, channel-to-field
mapping, conductor/electrode definitions and field instance. Uniform input
requires dimensioned `conductivityIntracellular` and
`conductivityExtracellular`. Field mode reads canonical
`ConductivityIntracellular` and `ConductivityExtracellular` using the existing
legacy handling. The monodomain effective tensor cannot uniquely identify that
physical pair and is not silently substituted for it.

Online projection copies voltage into `leadFieldObservationVm` and republishes
its source owner's sealed-wall policy under that private identity. The scheme
selector `leadFieldGradVm` is separate from the gradient/cache identity.
Neither source voltage, its boundaries nor a registered source gradient is
updated. A missing published policy is rejected.

The current metadata format is 2. Native sensitivities are checked against
physical snapshots and construction settings. Geometry fingerprints preserve
explicit myocardial restriction; full-conductor fingerprints additionally
check neighbours and physical patch names/types. Myocardial connectivity
fingerprints detect changed heart adjacency after heart-only export. Processor
connections use scalar hash chunks, preserving decomposition independence.

Heart-only import retains historical exterior provenance. It cannot recheck
an absent torso or its exterior boundaries. The fingerprints do not cover
all coupled-patch transformation options; AMI/nonconformal transfer and
rotated-cyclic configurations require separate validation. Moving/changing
meshes are rejected. Format-1 sensitivities must be regenerated.

Sensitivity filenames must be unique and cannot overwrite reserved source,
conductivity, metadata or helper fields. Validation happens before production
or synthetic construction and again at import.

## Tutorials and variants included in this PR

| Case | Location | Role |
|---|---|---|
| Monodomain manufactured observation | [monodomainLeadFieldECG](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/README.md) | Source voltage plus independent signed vector functional. |
| Bidomain manufactured observation | [bidomainLeadFieldECG](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/README.md) | Vm, gauge-corrected phiE and the independent functional. |
| Passive embedded-heart column | [leadFieldColumn](../tutorials/configurationChecks/base/leadFieldColumn/Allprepare) | Small reusable geometry with a nonzero construction/recording oracle. |
| Configuration variants | [configuration checks](../tutorials/configurationChecks/README.md) | Both source models, native/legacy conductivity fields, MPI, export, guards and recovery. |

The two manufactured cases are deliberately separate tutorials because they
have different source equations and mandatory field diagnostics. Their
optional pseudo-ECG is a runner variant, rather than another duplicated case:

```bash
./Allclean
./Allrun -pseudoECG yes
./Allclean
./Allrun fine -pseudoECG yes
./Allclean
./Allrun fine -dt .00005 -pseudoECG yes
```

`-N`, `-dt`, `-endTime` and `-pseudoECG` generate local includes without rewriting
tracked dictionaries. Regression runners exercise both observers in parallel.
`checkDualECG.py` checks common actual sample times, finite nonzero reference
cohorts, independent errors and subordinate reference-quadrature discrepancy.

Additional maintained post hooks cover input guards, heart-only export,
parallel construction, static forward recovery and rounded-folder restart.
The generic configuration harness supplies isolated run copies. Generated
meshes, fields, run logs and parameter includes are not PR inputs.

## What the verification exercises

The synthetic Vm is the existing FDA cosine-product reference. The added
sensitivity contains all three matching sine/cosine gradient components; its
exact volume integral is `-14*pi*kappa*A(t)/8`. Opposite kappa signs check
polarity and normalization. The verifier does not use production gradients or
cell sums to calculate that reference. Source voltage, gauge-corrected
extracellular potential and gradient errors are measured separately.

The mono fixture keeps the conductivity required by its fixed FDA forcing.
Its explicit observation pair Gi=Ge=2Gm has the intended monodomain reduction.
This corrects an earlier new-fixture mismatch that produced a voltage floor;
existing pseudo-ECG tutorial scientific defaults were not retuned.

The selected 20^3/40^3 pairs show LF slopes near 2.134, source Vm slopes near
2.022/1.998 and gradient slopes near 1.396. Fine-grid timestep perturbations
are small in the reported LF quantities. These remain two-level supporting
controls, not a settled asymptotic range or temporal-order result.

Concurrent pseudo-ECG uses three off-axis electrodes and independent
Gauss-Legendre reference orders 96/48. It uses the explicit raw prefactor one:
its amplitude is not compared directly with LF volts. Separate reference-
normalized cohort errors give LF 3.59% to 0.819% and pseudo 33.4% to 6.4%,
with slopes about 2.134/2.383. Different channel cohorts, kernels and
conductivity contracts prevent an intrinsic accuracy ranking.

The passive column independently gives Lx=50 m^-2 from unit current and area.
A frozen Vm=.01*x control gives 0.00475 V from both production recovery/torso
sampling and projection. It advances no ionic state or physical cardiac time.
It is a scalar orthogonal static comparison, not general torso equivalence.

Affine controls deliberately supply a stale registered gradient and check
unchanged voltage/extracellular snapshots. Rotated repeated-eigenvalue and
rank-one/two algebra controls exercise the SVD reconstruction. Serial/MPI,
field naming, dimensions, provenance, changed patch/connectivity, export and
restart checks supply separate software-contract evidence.

## Subsequent tutorials

A matched idealized heart/torso recording should be a variant of the existing
idealized-heart/bath workflow once the bath branch and electrode coordinates
are reconciled. It should share one geometry, prescribed source, conductor and
balanced leads between reciprocal and direct recordings. The current private
Niederer torso dataset is not bundled or treated as verified electrode input.

A nonorthogonal or rotated-anisotropic conductor MMS and a longer observation
mesh ladder would strengthen spatial evidence. They should extend the small
conductor fixtures before an anatomical accuracy claim. A rotated conformal
cyclic fixture would specifically test the hash-exchange path.

Eikonal/template support remains a separate extension. The current runtime
requires Vm from monodomain or bidomain; the existing eikonal pseudo-ECG test
does not verify that future reciprocal observer.

The final build/regression report records the tested revision and platform.
Cross-version/full-solids4foam CI and later application verification remain
separate from the local lightweight controls.

[The complete file inventory](LEAD_FIELD_ECG_FILES.md) lists every added or
extended path in the feature PR. [The final validation report](LEAD_FIELD_ECG_PR_VALIDATION.md)
records the tested source and remaining CI/application scope.
