# Lead-field ECG status and follow-ups

Status on 2026-10-13. Numerical implementation and maintained local controls
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
2. Add the new uniform-template eikonal observer and its selected-resolution
   result to Paper 1. Earlier monodomain/bidomain lead-field and dual-pseudo
   content is already integrated into Papers 1/2; the new eikonal result is
   not yet integrated. Paper 2 retains its bidomain/bath scope and matching
   methods/verification/results structure.
3. Refresh manuscript source/build provenance to the committed feature,
   archive the new raw traces/diagnostics immutably, and freeze the publication
   source and evidence release. Existing paper notes still include older
   uncommitted/superseded source descriptions.

The manuscripts are a separate working repository. Their other current edits
are preserved; they are not bundled in this code PR.

## Scientific/application follow-ups

These extend the verified scope; they are not prerequisites for the current
implementation PR.

- **Matched heart/torso recording:** extend the existing bath/idealized-heart
  workflow with the same source, physical conductor and electrode reference.
  Compare reciprocal recording with actual torso potential sampling through
  a full time series. Current torso recovery is a static scalar column control.
- **Anatomical Wilson leads:** verify E1–E6 and RA/LA/LL locations in the user's
  Niederer dataset, construct the corresponding reciprocal fields, and exercise
  them. Wilson weights are implemented; private coordinates are not certified
  inputs or bundled tutorial geometry.
- **Refinement for convergence-order claims:** monodomain/bidomain reciprocal
  evidence currently includes a 20³/40³ pair with observed order about 2.13.
  This is a supporting pair, not a settled asymptotic range. New eikonal
  reciprocal evidence is an accuracy check at 40³ (1.68% normalized RMS error).
  Add finer levels only when needed to support stronger paper claims, keeping
  source, sensitivity and sampling definitions fixed and recording field,
  gradient and functional errors separately. The 20³ eikonal constant-sensitivity
  integral is under-resolved and should not define an asymptotic slope.
- **Conductor reconstruction:** a rotated-anisotropic/nonorthogonal conductor
  case would extend the scalar orthogonal recovery evidence. Existing
  source/gradient shell controls do not automatically verify that conductor
  construction or the whole recording operator.
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
