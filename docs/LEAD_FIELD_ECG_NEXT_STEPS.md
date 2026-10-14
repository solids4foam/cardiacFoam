# Lead-field ECG status and follow-ups

Status on 2026-10-14. Numerical implementation and maintained local controls
are complete for monodomain, bidomain and spatially uniform-template eikonal
recording. The [validation report](LEAD_FIELD_ECG_PR_VALIDATION.md) identifies
the tested numerical source and platforms; the
[file inventory](LEAD_FIELD_ECG_FILES.md) lists the feature diff against PR #59.

## Completed scope

- Passive monodomain/bidomain lead-field observation with explicit physical
  Gi/Ge, static conductor preprocessing, Wilson channel definitions and
  geometry/topology/conductivity provenance checks.
- Simultaneous pseudo-ECG and reciprocal recording, each checked against its
  own reference; unchanged source fields and source-gradient ownership.
- Uniform-template eikonal recording through a shared template helper and
  gradient projection, including manufactured/endo/mid/epi unit/state controls.
- Maintained manufactured tutorials, static column recovery with actual torso
  sampling, MPI/import/export/input/restart controls, and implementation guides.
- Lightweight v2412/v2312 builds and selected regression sweeps. The full
  repository/full-solids4foam CI matrix belongs to the PR review stage.

## Release and paper integration

1. Review the stacked PR against `feat/mono-bi-eikonal-verification` (PR #59),
   inspect repository CI, and resolve applicable failures before merge.
2. Paper 1 now includes the uniform-template eikonal observer, N=10/20/40/80
   resolution controls and independent exact-source quadrature. Earlier
   monodomain/bidomain lead-field and dual-pseudo content is integrated into
   Papers 1/2. Paper 2 retains its bidomain/bath scope. The owner-requested
   serial N=80 extension is complete and independently checked.
3. Refresh manuscript source/build provenance to the committed feature,
   archive the new raw traces/diagnostics immutably, and freeze the publication
   source and evidence release. Existing paper notes still include older
   uncommitted/superseded source descriptions.

The manuscripts are a separate working repository. Their other current edits
are preserved; they are not bundled in this code PR.

## Scientific/application follow-ups

These extend the verified scope; they are not prerequisites for the current
implementation PR.

- **Matched heart/torso recording:** the new opt-in example now exercises
  a fresh 40-ms monodomain source, a short fresh heart-only bidomain source,
  and replay of an older coupled source. It records both online ECG models
  and independently recovers passive torso potentials on the same saved Vm.
  The measured monodomain forward/reciprocal cohort difference is 8.686%;
  selected tighter recovery controls do not explain it. Review the corrected
  operator/source-gradient relation and mesh sensitivity next; do not label
  the workflow's execution-contract pass an accuracy or equivalence result.
  A fresh full-window bidomain source run is still open. See the
  [application report](IDEALIZED_HEART_ECG_VERIFICATION.md).
- **Torso Wilson leads:** the six existing chest sites and declared torso
  equivalents for absent RA/LA/LL limbs now ship in the example. Native
  ownership, balanced referencing and independent Wilson reconstruction are
  exercised. The [matched-torso guide](MATCHED_HEART_TORSO_ECG.md) gives their
  coordinates. Its frame differs from the Niederer meshes; patient-specific
  registration and clinical lead validation remain separate.
- **Refinement for convergence-order claims:** monodomain/bidomain reciprocal
  evidence currently includes a 20³/40³ pair with observed order about 2.13.
  This is a supporting pair, not a settled asymptotic range. New eikonal
  reciprocal evidence includes N=10/20/40/80 resolution controls; N=10/20 are
  under-resolved, while N=40/80 have 1.679%/0.394% normalized RMS errors.
  The fine-pair reduction is 4.26; its iterative error has no tighter-criterion
  bound and a settled asymptotic range is not established.
  Add finer levels only when needed to support stronger paper claims, keeping
  source, sensitivity and sampling definitions fixed and recording field,
  gradient and functional errors separately. The 20³ eikonal constant-sensitivity
  integral is under-resolved and should not define an asymptotic slope.
- **Conductor reconstruction:** the idealized example now exposes the
  finite-mesh discrepancy on an anisotropic, nonorthogonal conductor with
  conductivity jumps. Its optional myocardial point fallback restores affine
  rank/exactness in 122 deficient face stencils; serial/six-rank checks pass.
  This does not establish general conductor accuracy or whole-operator
  equivalence. Existing source/gradient shell controls do not supersede that
  question. The original bath case's conormal interface model remains a
  separate control from the matched example's distance-weighted harmonic model.
- **Temporal evidence:** selected half-step controls already exist for the
  voltage-based observations. A formal temporal-order claim would need its own
  controlled study. Eikonal ECG sampling refinement concerns sampling/template
  interpolation rather than a propagation time step.
- **Regional eikonal templates:** specify and verify the spatial-weight term
  `sum_r U_r grad(w_r)` for smooth blending, and an interface interpretation for
  discontinuous regions. The reciprocal observer presently rejects this mode;
  legacy regional pseudo-ECG behavior remains available.

A complete discrete adjoint, anatomical/clinical equivalence and general
coupled-transform support are outside the current claims.
