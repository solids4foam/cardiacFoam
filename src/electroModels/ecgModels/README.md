# ecgModels

The ECG models, used by `ecgDomain`. They read the electrical state after each step and never feed back into the tissue.

## What's available

- `pseudoECG` (class `pseudoECGSolver`): a pseudo-ECG from the myocardium, using a dipole lead field (Gima and Rudy, 2002, Circulation Research, [doi:10.1161/01.res.0000016960.61087.86](https://doi.org/10.1161/01.res.0000016960.61087.86)).
- `torsoECG`: the potential at electrodes in the torso, the region outside the heart given by the top-level `torso` block of `electroProperties`. `torsoCoupling` selects how it is obtained:
  - `bath`: the global `phiE` of the heart-and-bath solve (`extracellularPotentialDomain`, bidomain only), sampled in the electrode cells;
  - `leadField`: the lead fields precomputed by `buildLeadFields` for a monodomain or bidomain heart, one-way. See [leadFieldECG](leadFieldECG/README.md).
- `eikonalECG`: a surrogate ECG for eikonal runs, built from activation times and action-potential templates. See [eikonalECG](eikonalECG/README.md).

The ECG domains are the entries of the top-level `ecgDomains` block of `electroProperties`; each entry accepts the keys of its `ecgSolver` only.

`pseudoECG` and `eikonalECG` share `infiniteMediumLeadField`: the lead field of an infinite homogeneous medium, scaled by `1/(4 pi sigmaExtracellular)`. Both require `sigmaExtracellular`, a positive scalar, in the ECG domain entry.

## Folders

```text
src/electroModels/ecgModels/
├── infiniteMediumLeadField/
├── pseudoECGSolver/
├── torsoECG/
├── leadFieldECG/
└── eikonalECG/
```
