# Immersed cylinder (static) smoke test

A deliberately tiny `pimpleHFDIBFoam` case whose job is to prove that the
immersed-boundary machinery still runs end to end. It is a smoke test, not a
validation case.

## What it is

A 2D channel, `x = [0, 24]`, `y = [0, 4.1]`, one cell thick in `z`, meshed by
`blockMesh` into **1250 cells**. A cylinder of diameter 1 centred at `(2, 2)` is
imposed as a static immersed body from `constant/triSurface/cylinder.stl`.
Inlet is a Poiseuille profile with `Umean = 1`; with `nu = 0.05` and `D = 1`
that is a nominal `Re = 20`.

It is derived from
`tutorials/pimpleHFDIBFoam/Immersed_Cyl_Static/Cyl_Immersed_New/ImmersedCyl-Static-Re20-Mesh1250`,
with three changes:

1. `endTime` reduced from 10 to 0.02, i.e. 20 steps of `deltaT = 0.001`.
2. The inlet `codedFixedValue` replaced by the identical profile written out as
   a `fixedValue nonuniform` list, so no run-time compilation is needed.
3. Initial conditions tracked as `0.orig/`, restored by `Allrun`.

## Running it

From a built repository:

```bash
cd tutorials/pimpleHFDIBFoam/immersedCylinderSmoke
./Allrun
```

As part of the regression suite:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression
```

`./Allclean` restores the case to its tracked state. Runtime is a few seconds
in serial; it needs no MPI, no plotting, no ParaView and no Python packages
beyond the standard library.

## What is asserted

`regression/checkSmoke.py` makes 16 checks, in three groups:

- **The run completed.** End time reached, log ends with `End`, no `nan`,
  `inf` or floating-point exception reported.
- **The immersed body exists and acts.** Exactly one active immersed body;
  `lambda` (solid volume fraction) finite and within `[0, 1]`; at least one
  cell more than half solid; total occupancy in a sane band; the immersed
  forcing field `f` finite, non-zero and bounded.
- **The solution has not diverged.** `U` and `p` finite and bounded, and the
  cumulative continuity error below `1e-4`.

Limits are wide bands around values measured on OpenFOAM v2412 and v2512. They
are chosen to catch a solver that crashes, diverges, produces `NaN` or quietly
stops forcing the body — not to pin down numbers.

## What it does *not* test

**No drag validation, and none is possible from this case.** The mesh grades
finest over `x = [10, 14]`, but the cylinder sits at `x = 2`, in the coarse
inlet region. It is resolved by roughly one cell across in `x` and six in `y`,
and only six cells carry any solid fraction at all. Any drag coefficient
computed here would be meaningless. This mesh/geometry mismatch is inherited
from the parent case and is present in every mesh level of that family
(1250 through 320000 cells), all of which place the cylinder at `(2, 2)`.

Also not covered: prescribed or force-driven body motion (this body is static,
`bodyOperation 0`), parallel decomposition, restart, and any of the DEM
contact, adhesion or virtual-mesh machinery.
