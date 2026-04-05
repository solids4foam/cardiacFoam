# verificationModels

Shared analytical manufactured-oracle definitions live here.

Current split:

- `monodomainVerification/manufacturedFDAMonodomainVerifier.C` owns the
  manufactured monodomain initialization and error-summary logic.
- `monodomainVerification/manufacturedFDAReference.H` provides the exact manufactured
  `Vm`, `u1`, `u2`, `u3`, and manufactured pseudo-ECG reference formulas.
- `bidomainVerification/manufacturedFDABidomainReference.H` extends the shared
  manufactured formulas with exact `phiE` and `phiI` fields for bidomain tests.
- `bidomainVerification/manufacturedFDABidomainVerifier.C` provides the
  bidomain manufactured error summary for `Vm`, `phiE`, `phiI`, `u1`, and `u2`.
- `ecgVerification/pseudoECGManufacturedVerifier.C` owns the manufactured
  pseudo-ECG verification lifecycle.
- Legacy compatibility wrappers in `electroModels` were removed; consumers
  should include `monodomainVerification/manufacturedFDAReference.H` directly.

The runtime-selected `tmanufacturedFDAPrePostProcessor` in
`modelPrePostProcessors/` is now only an adapter for the monodomain verifier.
