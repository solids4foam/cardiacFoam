# Bath domain: multi-mesh inventory

The bath domain introduced support for solving on a mesh that is distinct from
the main heart mesh.  Two strategies were implemented and selected at runtime:

| Strategy | Trigger | Class used |
|---|---|---|
| **Submesh** | `cellZone` key in dict | `fvMeshSubset` |
| **Region mesh** | `region` key in dict | owned `fvMesh` |

The files below are **not compiled** in this branch but remain on disk.  This
note records every class, member, and free function that carries multi-mesh
logic so the feature can be reinstated or ported cleanly.

---

## Classes

### `BathDomain`
`electroDomains/bathDomain/bathDomain.H` / `.C`

Implements three interfaces simultaneously:

- `electroDomainInterface` — plugs into the advance loop
- `bathCouplingEndpoint` — exposes fields to `heartBathInterfaceCoupler`
- `electroStateProvider` — lets the ECG solver read heart state

**Multi-mesh members**

| Member | Type | Purpose |
|---|---|---|
| `ownedMeshPtr_` | `autoPtr<fvMesh>` | Owned region mesh (non-null when strategy = region) |
| `meshSubsetPtr_` | `autoPtr<fvMeshSubset>` | Submesh carved from `supportMesh_` (non-null when strategy = submesh) |
| `supportMesh_` | `const fvMesh&` | Original heart mesh passed in at construction |
| `mesh_` | `const fvMesh&` | Reference to whichever mesh is active (resolved by `resolveBathMesh`) |
| `interfacePatchIndex_` | `label` | Patch index on the bath mesh that faces the heart |

**Multi-mesh methods**

| Method | Returns | Notes |
|---|---|---|
| `mesh()` | `const fvMesh&` | Active bath mesh (submesh or region) |
| `baseMesh()` | `const fvMesh&` | Always `supportMesh_` — the heart mesh |
| `subsetFaceMapPtr()` | `const labelUList*` | Non-null only for submesh; maps bath face indices → base mesh face indices (`fvMeshSubset::faceMap()`) |
| `subsetCellMapPtr()` | `const labelUList*` | Non-null only for submesh; maps bath cell indices → base mesh cell indices (`fvMeshSubset::cellMap()`) |
| `usesSubsetMesh()` | `bool` | True when `meshSubsetPtr_` is valid and has a sub-mesh |
| `interfacePatchIndex()` | `label` | Index of the heart-facing patch on the bath mesh, or -1 |

**Free functions** (file-local, `bathDomain.C`)

| Function | Purpose |
|---|---|
| `createBathMeshSubset(supportMesh, dict)` | Reads `cellZone` from dict, calls `fvMeshSubset::setCellSet` to carve the subset |
| `createBathRegionMesh(supportMesh, dict)` | Reads `region` from dict, constructs a new `fvMesh` for that region |
| `resolveBathMesh(supportMesh, ownedMeshPtr, meshSubsetPtr)` | Returns `ownedMeshPtr->mesh()`, `meshSubsetPtr->subMesh()`, or `supportMesh` depending on which pointer is valid |
| `bathPotentialPatchTypes(mesh, interfacePatchName)` | Builds the `wordList` of patch types for the `phiE` field, setting the interface patch to `fixedValue` |

---

### `BathECGSolver`
`ecgModels/bathECGSolver.H` / `.C`

Run-time selectable base class for solvers that act on a `BathDomain`.  No
multi-mesh logic of its own, but receives the active `BathDomain` (which may be
a submesh or region mesh) on each `solve()` call.

Concrete implementation: `bidomainBathECGSolver` in
`ecgModels/bidomainBathECGSolver/`.

---

### `bathCouplingEndpoint` (interface)
`electroCouplers/electroDomainCouplingEndpoints.H`

Abstract interface still present in the header (included by
`conductionSystemDomain.H`, `myocardiumDomainInterface.H`,
`electroDomainCoupler.H`).  **No compilation unit** — all methods are pure
virtual or have inline defaults.

| Method | Default | Purpose |
|---|---|---|
| `interfaceSourceField()` | pure virtual | Source term injected by coupler |
| `bathPotentialField()` | pure virtual | Solved extracellular potential on bath mesh |
| `interfacePatchIndex()` | returns `-1` | Patch index on bath mesh facing the heart |
| `baseMesh()` | pure virtual | Heart mesh from which bath was derived |
| `subsetFaceMapPtr()` | returns `nullptr` | Face map when submesh strategy is active |
| `subsetCellMapPtr()` | returns `nullptr` | Cell map when submesh strategy is active |

---

### `heartBathInterfaceCoupler`
`electroCouplers/heartBathInterfaceCoupler.H` / `.C`

Bridges `tissueCouplingEndpoint` (heart side) and `bathCouplingEndpoint` (bath
side).  Contains all the mapping arithmetic that transfers heart `phiE` values
onto the bath-mesh interface patch.

**Multi-mesh branching logic** (`preparePostPrimaryCoupling`, `.C`)

The coupler selects one of three transfer paths at runtime:

1. **Submesh path** — `subsetFaceMapPtr()` and `subsetCellMapPtr()` are
   non-null, and `baseMesh()` matches the heart mesh.  Walks `faceMap` to find
   the heart face that each bath interface face came from, then looks up the
   adjacent heart cell value.

2. **Region-mesh patch path** — `subsetFaceMapPtr()` is null,
   `interfacePatchIndex() >= 0`, and `baseMesh()` matches the heart mesh.
   Uses `treeDataFace` (nearest-face search) to match bath patch face centres
   to heart patch face centres and copies values.

3. **Same-mesh fallback** — bath mesh and heart mesh are identical;
   `phiE` primitive field is copied directly.

---

## Make file entries removed

The following lines were removed from `electroModels/Make/files` and
`electroModels/Make/options` in this branch.  Reinstate them when bath support
is added back.

**Make/files**
```
electroDomains/bathDomain/bathDomain.C
ecgModels/bathECGSolver.C
electroCouplers/heartBathInterfaceCoupler.C
```

**Make/options** (`EXE_INC`)
```
-IelectroDomains/bathDomain
```
