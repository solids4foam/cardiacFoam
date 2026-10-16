# Sealed-wall gradient probe

This test evaluates the `sealedWallLeastSquares` cell gradient of a prescribed
field on the monodomain `sphericalShell` meshes. It needs no ionic model and
no time step.

The field is either affine or the shell manufactured solution `F`
(`manufacturedShellMonodomainReference.H`). The probe publishes the
wall-gradient states of `psi` itself, as a solver's field owner does:
`endocardium` and `epicardium` sealed by default (`-sealed`), patches named
in `-value` as value patches, every other patch stock. A sealed face carries
the adjacent cell value, as a zero-gradient wall does; every other
non-coupled face carries the exact value.

Each run prints:

- `stockDiff`: the scheme with every patch stock against `leastSquares`;
- `nameDiff`: `grad(psi)` against `grad(psi.component(0))` evaluated on a
  copy of the field under another name;
- `wallRMS`, `wallMax`: gradient error in the cells with a sealed or value
  face, relative to the largest exact gradient;
- `interiorRMS`, `interiorMax`: the same in the other cells.

Run from any directory after setting up OpenFOAM:

```bash
applications/utilities/wallGradientProbe/run.sh hex 8 16 32
```

```bash
applications/utilities/wallGradientProbe/run.sh tet 6 10 14
```

`run.sh` meshes each resolution (`blockMesh` cubed sphere or `gmsh`
tetrahedra), runs the probe serial and on 4 scotch ranks, and fails unless:

- `stockDiff` and `nameDiff` are exactly zero;
- the affine field's gradient error is below `1e-10` in every cell;
- the serial and parallel shell errors agree to `1e-10` relative.

It then prints the shell errors with the observed order between levels.
The wall-cell gradient is first order on fine enough meshes; the interior
is the stock `leastSquares` gradient.
