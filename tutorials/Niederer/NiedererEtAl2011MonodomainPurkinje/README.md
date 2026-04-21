# NiedererEtAl2011Purkinje tutorial architecture

This tutorial keeps the Niederer slab setup and adds a small 1D Purkinje
network coupled into the 3D monodomain tissue.

- Myocardium solver: `monodomainSolver`
- Conduction solver: `monodomain1DSolver`
- PVJ coupler: `reactionDiffusionPvjCoupler`
- Goal: document a minimal Purkinje-enabled slab case without changing the core
  Niederer probe layout

## Folder structure

```text
tutorials/NiedererEtAl2011Purkinje/
├── constant/
│   ├── electroProperties
│   └── physicsProperties
├── system/
│   ├── Niedererlines
│   ├── Niedererpoints
│   ├── blockMeshDict
│   ├── controlDict
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── setupNiedererEtAl2011Purkinje/
│   └── README.md
├── Allrun
├── Allclean
└── README.md
```

## Key dictionary additions

`constant/electroProperties` extends the base Niederer slab with:

- `conductionNetworkDomains.purkinjeNetwork` for a small 4-node network
- `purkinjeNetworkModelCoeffs.conductionSystemSolver monodomain1DSolver`
- `domainCouplings.couplingA` using `reactionDiffusionPvjCoupler`

The tutorial keeps `couplingMode unidirectional` so the standard staggered
advance scheme remains valid for a simple explanatory case.

## Run modes

```bash
./Allrun
./Allrun parallel
```

Use the existing Niederer point and line probes to inspect how the Purkinje
terminals alter the activation pattern near the PVJ locations.
