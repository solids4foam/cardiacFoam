# configurations - bathBidomain/obliqueSlab

The interface configurations 0, A, AI, ATI and ATBI (see the case README) at
fixed resolution, on hexahedra and tetrahedra.

## Execution

In a copy of the case, for each configuration:

```bash
./Allrun parallel -N 40 -config <config>
./Allrun parallel -N 40 -config <config> -elliptic -nOuter 10 -nNonOrth 10
./Allrun parallel -mesh tet -N 20 -config <config>
```

Compare `Vm`, `phiE_heart`, `phiE_bath`, `interfaceFlux` and the wall rows
across configurations; the elliptic orders by configuration are in the case
README's expected outcomes.

## Outputs

Run copies, meshes and results stay outside the repository (`results/` is
ignored).
