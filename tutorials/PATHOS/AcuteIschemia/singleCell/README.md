# singleCell

Single integration-point electrophysiology case for the AcuteIschemia pathology.
Uses the same `ionicConstantOverrides` as the parent tissue case to validate the
modified action potential shape before a full 3D run.

- Electro model: `singleCellSolver`
- Ionic model: `TNNP`
- Ionic overrides: hyperkalemia + INa/ICaL/IK1/INaK modifications (from `ionicConstantOverrides`)

## Folder structure

```text
tutorials/PATHOS/AcuteIschemia/singleCell/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── studies/
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
