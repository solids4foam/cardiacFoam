# purkinjeNiedererEtAl2011 tutorial architecture

This tutorial keeps the Niederer slab setup and adds a small 1D Purkinje
network coupled into the 3D monodomain tissue.

- Myocardium solver: `monodomainSolver`
- Conduction solver: `monodomain1DSolver`
- PVJ coupler: `reactionDiffusionPvjCoupler`

## Folder structure

```text
tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/
├── constant/
│   ├── electroProperties
│   ├── physicsProperties
│   └── purkinjeGraph
├── system/
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── Allrun
├── Allclean
├── purkinjeSlab.reference
├── regressionTest.sh
└── README.md
```

## Key dictionary additions

`constant/electroProperties` extends the base Niederer slab with:

- `conductionNetworkDomains.purkinjeNetwork` for a small 4-node network
- `purkinjeGraphModelCoeffs.conductionSystemSolver monodomain1DSolver`
- `constant/purkinjeGraph` for edges, points, PVJ nodes, and PVJ locations
- `domainCouplings.couplingA` using `reactionDiffusionPvjCoupler`

The tutorial keeps `couplingMode unidirectional` so the standard staggered
advance scheme remains valid for a simple explanatory case.

## Run modes

```bash
./Allrun
./Allrun parallel
```
