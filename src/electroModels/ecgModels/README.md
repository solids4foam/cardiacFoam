# ecgModels

The ECG models, used by `ecgDomain`. They read the electrical state after each step and never feed back into the tissue.

## What's available

- `pseudoECG` (class `pseudoECGSolver`): a pseudo-ECG from the myocardium, using a dipole lead field (Gima and Rudy, 2002, Circulation Research, [doi:10.1161/01.res.0000016960.61087.86](https://doi.org/10.1161/01.res.0000016960.61087.86)).
- `torsoECG`: samples the global `phiE` of a heart-and-bath solve at the electrodes.
- `eikonalECG`: a surrogate ECG for eikonal runs, built from activation times and action-potential templates. See [eikonalECG](eikonalECG/README.md).

## Folders

```text
src/electroModels/ecgModels/
├── pseudoECGSolver/
├── torsoECG/
└── eikonalECG/
```
