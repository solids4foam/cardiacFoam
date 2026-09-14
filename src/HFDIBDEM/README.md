# HFDIBDEM library: provenance and attribution

This directory contains a **modified derivative** of the third-party
`openHFDIB-DEM` library. It is not original cardiacFoam code. This note records
where the code came from, what was changed, and which points are still
unconfirmed.

This is a factual provenance record, not legal advice.

## Upstream project

| Item | Value |
| --- | --- |
| Project | openHFDIB-DEM |
| Primary repository | <https://github.com/techMathGroup/openHFDIB-DEM> |
| Earlier author fork | <https://github.com/MartinIsoz/openHFDIB-DEM> |
| Repository licence | GNU GPL v3 (`LICENSE` in the upstream repository) |
| Reference paper | Studeník, O.; Isoz, M.; Kotouč Šourek, M.; Kočí, P. *OpenHFDIB-DEM: An extension to OpenFOAM for CFD-DEM simulations with arbitrary particle shapes.* SoftwareX **27** (2024) 101871. DOI: [10.1016/j.softx.2024.101871](https://doi.org/10.1016/j.softx.2024.101871) |
| Method papers | Municchi & Radl (2016); Isoz & Kotouč Šourek (2020), as cited in the file headers |

Upstream contributors named in the preserved file headers: Federico Municchi
(2016), Martin Isoz (2019–), Martin Kotouč Šourek (2019–), Ondřej Studeník
(2020–).

If you use this library, please cite the SoftwareX paper above.

## Licence status

- The upstream repository ships a **GPL v3** `LICENSE` file and its README
  states GPL v3.
- The upstream **per-file headers** instead say
  *"openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE
  (LGPL)"*.
- **This discrepancy exists upstream and was not introduced by cardiacFoam.**
  It was verified against the upstream repository at the revisions listed
  below.
- The file headers here are reproduced **verbatim** from upstream. They have
  not been rewritten, replaced or relicensed.
- cardiacFoam itself is GPL v3 (root `LICENSE`), so redistributing this code
  under GPL v3 is consistent with the upstream repository licence under either
  reading of the discrepancy.

Anyone needing certainty about the LGPL-vs-GPL question should raise it with
the upstream maintainers rather than editing these headers locally.

Because these headers are preserved verbatim, this directory is excluded from
the repository's file-header check; see the exclusion list and its rationale in
`.github/workflows/checkFileHeader.yml`. Do not run
`applications/scripts/cardiacFoamChangeCopyright.sh` on files here.

## Import path

The code did not come directly from upstream. It arrived through an
intermediate UCD repository:

```text
techMathGroup/openHFDIB-DEM          (upstream, GPL v3 repo / LGPL headers)
  └─ xenosim-erc/immersedBoundaryRigidMotion
       │  port to OpenFOAM v2412/v2512, removal of DEM particle,
       │  contact-model and virtual-mesh machinery, Clang fixes
       └─ cardiacFoam  src/HFDIBDEM  (this directory)
            valve-specific prescribed motion
```

Path mapping from upstream to here:

| Upstream path | Path here |
| --- | --- |
| `src/HFDIBDEM/**` | `src/HFDIBDEM/**` |
| `applications/solvers/incompressible/pimpleHFDIBFoam/**` | `applications/solvers/pimpleHFDIBFoam/**` |
| `pimpleHFDIBFoam/correctPhi.H` | `applications/solvers/pimpleHFDIBFoam/correctPhi.solver.H` |

## Known modifications

Relative to the intermediate repository, the following files carry
cardiacFoam-specific changes, mainly to support prescribed valve motion driven
by a per-vertex geodesic coordinate field:

`addModels/addModel.C`, `addModels/initializeAddModels.H`,
`geomModels/stlBased/stlBased.C`, `geomModels/stlBased/stlBased.H`,
`ibInterpolation/leastSquaresInt/leastSquaresInt.C`,
`ibInterpolation/leastSquaresInt/leastSquaresIntInfo.C`,
`ibInterpolation/lineInt/lineInt.C`, `ibInterpolation/lineInt/lineInt.H`,
`ibInterpolation/lineInt/lineIntInfo.C`,
`ibInterpolation/lineInt/lineIntInfo.H`, `immersedBody.C`, `immersedBody.H`,
`openHFDIBDEM.C`, `openHFDIBDEM.H`.

Relative to upstream, the intermediate repository additionally removed the DEM
contact models, particle add-models and virtual-mesh subsystem, and ported the
remainder to modern OpenFOAM.

`applications/solvers/pimpleHFDIBFoam` is a modified `pimpleFoam`. Its
`pimpleHFDIBFoam.C` now carries cardiacFoam's GPL v3 licence block with the
upstream `Copyright (C)` lines retained (OpenFOAM Foundation 2011–2017,
OpenCFD Ltd 2019) and the note "Ported by Philip Cardiff". Its `.H` include
fragments are left headerless, as upstream OpenFOAM ships them.

## Outstanding confirmations

These are recorded as unknown rather than guessed:

1. The exact upstream revision the original import was taken from. All
   unmodified files match upstream exactly, but 29 upstream revisions are
   equally consistent with them. The evidence is compatible with the `Ver2.6`
   branch and rules out anything older than 2024-05-01.
2. Whether the LGPL header text or the GPL v3 `LICENSE` is upstream's intended
   licence.
3. Attribution for the small number of files in this directory that carry no
   licence or author header (`declareExternVars.H`, `defineExternVars.H`,
   `initializeIB.H`, `addModels/initializeAddModels.H`).
