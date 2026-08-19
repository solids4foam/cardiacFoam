# singleCell

Single integration-point electrophysiology case for the BrugadaSyndrome pathology.
Uses the same `ionicConstantOverrides` as the parent tissue case to validate the
modified action potential shape before a full 3D run.

- Electro model: `singleCellSolver`
- Ionic model: `TNNP`
- Ionic overrides: INa −80% (global), Ito ×6 + ICaL −50% (epicardialCells)

## Folder structure

```text
tutorials/PATHOS/BrugadaSyndrome/singleCell/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── setup/
├── singleCell.reference
├── regressionTest.sh
├── Allrun
└── Allclean
```

## Run

```bash
./Allrun
./regressionTest.sh
```
