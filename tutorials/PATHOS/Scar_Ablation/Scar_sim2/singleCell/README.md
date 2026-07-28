# singleCell

Single integration-point electrophysiology case for the Scar_Ablation pathology.
Uses the same `ionicConstantOverrides` as the parent tissue case to validate the
action potential shape in scar and border zone configurations.

- Electro model: `singleCellSolver`
- Ionic model: `TNNP`

## Folder structure

```text
tutorials/PATHOS/Scar_Ablation/Scar_sim2/singleCell/
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
