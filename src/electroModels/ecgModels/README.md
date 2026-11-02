# ecgModels

The ECG models, used by `ecgDomain`. They read the electrical state after each step and never feed back into the tissue.

## What's available

- `pseudoECG` (class `pseudoECGSolver`): a pseudo-ECG from the myocardium, using a dipole lead field (Gima and Rudy, 2002, Circulation Research, [doi:10.1161/01.res.0000016960.61087.86](https://doi.org/10.1161/01.res.0000016960.61087.86)).
- `eikonalECG`: a surrogate ECG for eikonal runs, built from activation times and action-potential templates. See [eikonalECG](eikonalECG/README.md).
- `torsoECG`: the potential at electrodes in the torso, the region outside the heart. With a bath it samples the potential the bath solves with a bidomain heart every step; with lead fields it contracts lead fields precomputed for the torso with the heart's source, one-way, for a monodomain or bidomain heart. See [leadFieldECG](leadFieldECG/README.md).

`pseudoECG` and `eikonalECG` share `infiniteMediumLeadField`, the lead field of an infinite homogeneous medium; they need no torso mesh.

## Folders

```text
src/electroModels/ecgModels/
├── infiniteMediumLeadField/
├── pseudoECGSolver/
├── torsoECG/
├── leadFieldECG/
└── eikonalECG/
```
