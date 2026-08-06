# The `physicsModel` One Definition Rule (ODR) Violation

## The Problem (Bus Error 10)

During the initialization of `cardiacFoam`, the `electroMechanicalNiedererEtAl2011` tutorial (and any other solver) would crash immediately with a `Bus error: 10` after the mesh was created.

The crash occurred inside `electroMechanicalModel.C` when trying to access the solid mesh through `solid_->mesh()`.

## Root Cause

The root cause was an inconsistent class layout due to a subtle **C++ One Definition Rule (ODR) violation** across two different modules:

1. **`solids4foam`** defines a base class `physicsModel` in its own library (`libsolids4FoamModels.dylib`).
2. **`cardiacFoam`** overrides this header by providing a duplicate version of `physicsModel.H` inside `modules/physicsModel/src/solids4FoamModels/lnInclude/`.
3. In `cardiacFoam`'s version, two extra pointers were added to support FSI simulations:

   ```cpp
   autoPtr<dynamicFvMesh> fluidMeshPtr_;
   autoPtr<dynamicFvMesh> solidMeshPtr_;
   ```

   This made the `physicsModel` class **16 bytes larger** than the original version in `solids4foam`.

### Why the build system allowed this

Even though `solids4foam` is compiled inside the `cardiacFoam` folder, the OpenFOAM build system (`wmake`) treats `solids4foam` as an independent, self-contained library.

- When `solids4foam` compiles, its `Make/options` file only looks inside its own folders. It uses its internal original version of `physicsModel.H` and allocates `solidModel` objects **without** the two pointers.
- When `cardiacFoam` compiles, its `Make/options` explicitly includes the modified `modules/physicsModel` folder. The compiler assumes `solidModel` objects **have** the two pointers and are 16 bytes larger.

### The Crash

At runtime, `cardiacFoam` asks `solids4foam` to allocate a `solidModel`. `solids4foam` allocates the smaller object. When `cardiacFoam` tries to access `solid_->mesh()`, it calculates the memory offset based on the *larger* header. Because of this 16-byte mismatch, `cardiacFoam` reads past the actual variables and hits garbage memory, resulting in the `Bus error`.

## The Fix

Since `fluidMeshPtr_` and `solidMeshPtr_` were only initialized but completely unused throughout the rest of the `cardiacFoam` codebase, they were **removed** from the modified `physicsModel.H` and `physicsModel.C` files.

This successfully restored the ABI memory compatibility between `cardiacFoam` and the `solids4foam` shared libraries.

## Future FSI Implementation

To safely add those pointers for FSI without crashing, you cannot modify the base class that `solids4foam` relies on (unless you force `solids4foam`'s `Make/options` to use `cardiacFoam`'s headers, which is hard to maintain).

Instead, the standard C++ approach is to leave the shared `physicsModel` base class alone, and create a **brand new derived class** (e.g., `fsiPhysicsModel : public physicsModel`) strictly inside `cardiacFoam` where the extra pointers can be safely added.
