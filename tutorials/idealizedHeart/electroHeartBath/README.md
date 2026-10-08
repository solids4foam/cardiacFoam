# Bath bidomain on the idealized heart

The electroHeart monodomain case as a bidomain myocardium embedded in an
idealized torso: one conformal heart+torso mesh, a global extracellular
potential `phiE` solved through the heart and the torso, and the precordial
leads sampled in the torso.

## Stack

| part | setting |
|---|---|
| myocardium | `bidomainSolver` on cellZone `myocardium`, `solutionAlgorithm implicit` |
| torso | `bathPotentialDomain`, zones `blood lungs chestWall torso`, `distanceWeightedHarmonic` interface conductivity, `bathPredictorCorrector yes` |
| ionic | `BuenoOrovio`, endo/M/epi `namedRegions` on `t`, apicobasal `tauSi` gradient |
| conduction system | human `purkinjeGraph`, `monodomain1DSolver`, `reactionDiffusionPvjCoupler` |
| pre-pacing | `constant/prePacingProperties`, as in electroHeart |
| ECG | `torsoECG` (phiE in the torso) and `pseudoECG` (infinite medium, sigma 0.2) at the same V1-V6 |

## Mesh and fields

`Allrun` copies everything it reads from `../mesh/heartTorso`, the
heart+torso counterpart of `../mesh`, stored with LFS like the rest of
`../mesh`:

| file | content |
|---|---|
| `constant/polyMesh` | 299,482 tetrahedra; one patch `torsoSurface`; cellZones `myocardium blood lungs chestWall torso` |
| `0/fiber sheet sheetNormal tm tv apicobasal` | the shared heart's fields, bit-identical; 0 in the torso |
| `0/t` | `1 - tm`, 0 at the endocardium and 1 at the epicardium |
| `0/Conductivity`, `0/ConductivityIntracellular`, `0/ConductivityExtracellular` | cardiacCore `setCardiacConductivity` output; 0 in the torso |
| `0/bodyAndOrgansConductivity` | torso zone conductivities |
| `constant/torsoGeometry` | torso box and organ shapes in global coordinates |

The heart cells are cells 0..123,616 in the shared heart's order, so the
shared `../mesh/constant/purkinjeGraph` applies unchanged.

`meshGeneration/Allmesh` writes `../mesh/heartTorso` and is run only to
change the torso, the organs or the conductivities:

1. `heartMesh/` receives the shared heart (mesh and anatomy fields).
2. `build_torso_mesh.py` reads the heart boundary triangles (`EPI`, `BASE`,
   `ENDO_LV`, `ENDO_RV`) and gives them to gmsh as fixed discrete surfaces;
   gmsh tetrahedralises the torso box minus the heart (2 mm at the heart
   grading to 12 mm at the torso surface).
3. `gmshToFoam`, `mergeMeshes` (heart cells first), `stitchMesh` (four
   `perfect` patch pairs, `system/stitchMeshDict`), `createPatch` (drops the
   emptied patches).
4. `topoSet` builds the zones; `setTorsoOrganConductivityField` builds
   `bodyAndOrgansConductivity`.
5. `mapFieldsPar` (`cellVolumeWeight`) fills the heart cells of
   zero-initialised fields from `heartMesh/0`; `setExprFields` writes `t`;
   `setCardiacConductivity` writes the three conductivity tensors.

It needs the gmsh Python module and cardiacCore's `setCardiacConductivity`
on `PATH`.

## Torso geometry

The torso is a box in a frame rotated from the mesh axes by
`TORSO_EULER_DEG` (`meshGeneration/build_torso_mesh.py`): `e1` caudal,
`-e2` patient left, `-e3` anterior; the heart long axis (base to apex)
points caudal, left and anterior. The anterior face passes through the most
anterior of V1-V3, the left face through the most lateral of V4-V6. The
electrodes keep electroHeart's coordinates; V2 and V5 lie on the surface,
the others at most 8.5 mm beneath it, and the heart keeps at least 15 mm
from the wall. The box is 0.20 m (anterior-posterior) x 0.26 m (left-right)
x 0.32 m (cranio-caudal).

The organs are not meshed. They are analytic shapes in the torso frame; the
zones follow the tetrahedra. Only the heart-torso interface is conformal:
faces between two torso cells take a linear average of the cell
conductivities, so organ boundaries need no conformal surface.

| zone | shape | cells | sigma [S/m] |
|---|---|---|---|
| `blood` | torso cells apical of the base plane connected to an LV or RV cavity point (`regionToCell`) | 31,192 | 0.7 |
| `lungs` | two cylinders along `e1`, posterior-left and right of the heart | 9,680 | 0.0389 |
| `chestWall` | the outer 10 mm of the box: skin, fat, muscle and ribs as one layer | 17,027 | 0.05 |
| `torso` | every other torso cell | 117,966 | 0.22 |

Precedence is chestWall, blood, lungs; torso takes the remainder, so the
zones partition the torso cells as `bathPotentialDomain` requires. Blood,
lung and torso values are the nominal forward-ECG torso values tabulated in
arXiv:2407.17146 (Table 2); the chest wall is the series mean of fat
(0.037) and skeletal muscle (~0.2) lowered for skin (0.01) and bone
(0.006). `phiE` is 0 at `phiERefPoint`, 20 mm inside the
posterior-caudal-right corner.

## Bidomain conductivities

Johnston, "Six conductivity values to use in the bidomain model of cardiac
tissue", IEEE TBME 63(7), 2016, nominal set, in the fibre / sheet / normal
frame:

| | fibre | sheet | normal |
|---|---|---|---|
| intracellular [S/m] | 0.24 | 0.035 | 0.008 |
| extracellular [S/m] | 0.24 | 0.20 | 0.11 |
| monodomain-equivalent `gi ge/(gi + ge)` | 0.12 | 0.0298 | 0.0075 |

The electroHeart `Conductivity` is (0.5, 0.1, 0.1), so conduction here is
slower than in electroHeart.

## Usage

```bash
./Allrun                 # 40 ms (QRS), controlDict.qrs
./Allrun fullBeat        # 400 ms (QRS and T wave)
./Allrun parallel        # decomposePar / cardiacFoam / reconstructPar
./Allclean
meshGeneration/Allmesh   # rewrite ../mesh/heartTorso
meshGeneration/Allclean  # remove the generation intermediates
```
