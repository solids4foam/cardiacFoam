# checkMeshGeometry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a standalone OpenFOAM utility that reads `constant/polyMesh`, detects if points are not in SI meters by order of magnitude, warns the user, and auto-scales + rewrites the mesh if needed.

**Architecture:** Single `.C` file utility following the existing pattern in `applications/utilities/`. Reads a `polyMesh`, computes bounding box, applies threshold-based unit detection, scales in-place via `movePoints`, and rewrites with `mesh.write()`.

**Tech Stack:** OpenFOAM (polyMesh, boundBox, pointField), wmake

---

### Task 1: Create Make/files

**Files:**

- Create: `applications/utilities/checkMeshGeometry/Make/files`

**Step 1: Create the file**

```
checkMeshGeometry.C

EXE = $(FOAM_USER_APPBIN)/checkMeshGeometry
```

**Step 2: Verify it exists**

```bash
cat applications/utilities/checkMeshGeometry/Make/files
```

Expected: shows the two lines above.

---

### Task 2: Create Make/options

**Files:**

- Create: `applications/utilities/checkMeshGeometry/Make/options`

**Step 1: Create the file**

Only `finiteVolume` and `meshTools` are needed — `polyMesh` and `boundBox` are covered by these.

```
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

EXE_LIBS = \
    -lfiniteVolume \
    -lmeshTools
```

**Step 2: Verify it exists**

```bash
cat applications/utilities/checkMeshGeometry/Make/options
```

---

### Task 3: Implement checkMeshGeometry.C

**Files:**

- Create: `applications/utilities/checkMeshGeometry/checkMeshGeometry.C`

**Step 1: Write the implementation**

Unit detection logic uses explicit thresholds (more robust than pure log10 formula):

| maxDim range    | Assumed unit | scaleFactor |
|-----------------|--------------|-------------|
| < 1.0           | m (correct)  | 1 (no-op)   |
| 1.0 – 999       | mm           | 1e-3        |
| 1000 – 999 999  | μm           | 1e-6        |

```cpp
/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    checkMeshGeometry

Description
    Reads constant/polyMesh and checks whether the mesh points are in SI
    meters. If not, warns the user and rescales the mesh to meters based
    on the detected order of magnitude (mm -> 1e-3, um -> 1e-6).

Usage
    checkMeshGeometry

Author
    cardiacFoam developers.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "boundBox.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Check mesh geometry units and auto-scale to SI meters if needed."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== checkMeshGeometry ==========\n" << endl;

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );

    const boundBox bb(mesh.points(), false);
    const scalar maxDim = cmptMax(bb.span());

    Info<< "Bounding box : " << bb << nl
        << "Max dimension: " << maxDim << " m" << nl << endl;

    // Threshold-based unit detection
    scalar scaleFactor = 1.0;
    word detectedUnit = "m";

    if (maxDim >= 1.0 && maxDim < 1000.0)
    {
        scaleFactor = 1e-3;
        detectedUnit = "mm";
    }
    else if (maxDim >= 1000.0 && maxDim < 1e6)
    {
        scaleFactor = 1e-6;
        detectedUnit = "um";
    }

    if (mag(scaleFactor - 1.0) < SMALL)
    {
        Info<< "Mesh appears to be in meters. No scaling applied." << nl
            << endl;
    }
    else
    {
        WarningInFunction
            << "Max dimension = " << maxDim
            << " suggests mesh is in " << detectedUnit
            << ", not meters.\n"
            << "  Applying scale factor: " << scaleFactor << "\n"
            << "  Rewriting constant/polyMesh ..." << nl
            << endl;

        pointField newPoints(mesh.points());
        newPoints *= scaleFactor;
        mesh.movePoints(newPoints);
        mesh.write();

        const boundBox bbScaled(mesh.points(), false);
        Info<< "Scaled bounding box : " << bbScaled << nl
            << "Scaled max dimension: " << cmptMax(bbScaled.span()) << " m"
            << nl << endl;
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
```

---

### Task 4: Build the utility

**Step 1: Build**

```bash
cd applications/utilities/checkMeshGeometry
wmake
```

Expected: no errors, binary appears at `$FOAM_USER_APPBIN/checkMeshGeometry`.

**Step 2: Verify binary exists**

```bash
which checkMeshGeometry
```

Expected: prints a path ending in `checkMeshGeometry`.

**Step 3: Commit**

```bash
git add applications/utilities/checkMeshGeometry/
git commit -m "feat(checkMeshGeometry): add standalone mesh unit detection utility"
```

---

### Task 5: Smoke-test — mesh already in meters

Use the NiedererEtAl2012 tutorial (blockMesh with `scale 0.001` → mesh in meters, maxDim ≈ 0.02 m).

**Step 1: Generate mesh**

```bash
cd tutorials/NiedererEtAl2012
blockMesh
```

**Step 2: Run checkMeshGeometry**

```bash
checkMeshGeometry
```

Expected output (no scaling):

```
========== checkMeshGeometry ==========

Bounding box : (0 0 0) (0.02 0.003 0.007)
Max dimension: 0.02 m

Mesh appears to be in meters. No scaling applied.

End
```

---

### Task 6: Smoke-test — mesh in mm (wrong units)

**Step 1: Scale mesh to mm using OpenFOAM's transformPoints**

```bash
cd tutorials/NiedererEtAl2012
transformPoints -scale '(1e3 1e3 1e3)'
```

This simulates a mesh imported in mm (maxDim becomes 20.0).

**Step 2: Run checkMeshGeometry**

```bash
checkMeshGeometry
```

Expected output (warns + scales):

```
========== checkMeshGeometry ==========

Bounding box : (0 0 0) (20 3 7)
Max dimension: 20 m

--> FOAM Warning :
    From function int main(int, char**)
    ...
    Max dimension = 20 suggests mesh is in mm, not meters.
      Applying scale factor: 0.001
      Rewriting constant/polyMesh ...

Scaled bounding box : (0 0 0) (0.02 0.003 0.007)
Scaled max dimension: 0.02 m

End
```

**Step 3: Re-run checkMeshGeometry to confirm idempotency**

```bash
checkMeshGeometry
```

Expected: "Mesh appears to be in meters. No scaling applied."

**Step 4: Commit test results (no file changes expected — just verify)**

```bash
git checkout tutorials/NiedererEtAl2012/constant/polyMesh/points
```

Restore the correct points for the tutorial case.

---

### Task 7: Update Allwmake comment in docs (optional)

`applications/Allwmake` already uses `wmake all utilities` which picks up all subdirectories automatically — no change needed.

Verify:

```bash
cd applications
./Allwmake 2>&1 | grep checkMeshGeometry
```

Expected: shows wmake processing `checkMeshGeometry`.

---

### Task 8: Update MEMORY.md

Add a short entry for this utility:

```markdown
## checkMeshGeometry utility
- `applications/utilities/checkMeshGeometry/checkMeshGeometry.C`
- Reads `constant/polyMesh`, detects units by threshold on maxDim
- maxDim < 1.0 → m (no-op), 1–999 → mm (×1e-3), 1000–1e6 → μm (×1e-6)
- Call after `newVtkUnstructuredToFoam` in Allrun scripts
```

**Commit:**

```bash
git add .claude/
git commit -m "docs(memory): add checkMeshGeometry utility entry"
```
