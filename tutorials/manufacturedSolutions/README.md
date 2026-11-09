# Manufactured solutions

Verification cases that run a solver against a known (manufactured or exact)
solution and check the error against the reference solution. The cases are
grouped by the equation they solve, then by the geometry they exercise.

```text
manufacturedSolutions/
├── generators/       the definition of each manufactured solution
├── monodomain/       cube, sphericalShell, coupled1D3DCube, coupled1D3DShell
├── eikonal/          cube, sphericalShell, coupled1D3DSphere
├── bidomain/         cube, sphericalShell
├── bathBidomain/     slab
└── electromechanics/ totalLagrangianEM
```

## Cases

| Equation | Case | Geometry |
| --- | --- | --- |
| [monodomain](monodomain/README.md) | [`cube`](monodomain/cube/README.md) | unit cube, hexahedra or tetrahedra, pseudo-ECG |
| | [`sphericalShell`](monodomain/sphericalShell/README.md) | spherical shell with sealed curved walls |
| | [`coupled1D3DCube`](monodomain/coupled1D3DCube/README.md) | Purkinje graph joined to the unit cube |
| | [`coupled1D3DShell`](monodomain/coupled1D3DShell/README.md) | Purkinje graph joined to the shell, tetrahedra |
| [eikonal](eikonal/README.md) | [`cube`](eikonal/cube/README.md) | unit cube, hexahedra or tetrahedra, eikonal ECG |
| | [`sphericalShell`](eikonal/sphericalShell/README.md) | spherical shell with sealed curved walls |
| | [`coupled1D3DSphere`](eikonal/coupled1D3DSphere/README.md) | graph into a sealed tissue sphere, tetrahedra |
| [bidomain](bidomain/README.md) | [`cube`](bidomain/cube/README.md) | unit cube, hexahedra or tetrahedra |
| | [`sphericalShell`](bidomain/sphericalShell/README.md) | spherical shell with sealed curved walls |
| [bathBidomain](bathBidomain/README.md) | [`slab`](bathBidomain/slab/README.md) | myocardium slab between two bath slabs |
| [electromechanics](electromechanics/README.md) | [`totalLagrangianEM`](electromechanics/totalLagrangianEM/README.md) | unit cube, needs solids4foam |

The Niederer slab benchmark is in
[`../benchmarks`](../benchmarks/NiedererEtAl2011verification/README.md).

## Definitions of the manufactured solutions

[`generators/`](generators/README.md) holds the definition of each monodomain
and eikonal manufactured solution: the sympy generators, the NumPy reference
modules, and the two construction documents
([monodomain_mms.pdf](generators/docs/monodomain/monodomain_mms.pdf),
[eikonal_mms.pdf](generators/docs/eikonal/eikonal_mms.pdf)). The C++ reference
headers in `src/verificationModels/` are generated from the generators, in
whole or in marked regions. The case READMEs describe how to run a case and
point to these definitions.

## Run a case

Each case has its own `Allrun` (`./Allrun -h` lists its options), `Allclean`
and `regression/regressionTest.sh`. There is no suite-level `Allrun`.

```bash
cd tutorials/manufacturedSolutions/bidomain/cube
./Allrun -h                       # options of this case
./Allrun                          # mesh, run, post-process
./Allrun -mesh tet                # the same case on gmsh tetrahedra
```

`tutorials/Alltest-regression` runs the regressions; `--only` selects cases by
path:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'manufacturedSolutions/bidomain'
```

See [`../README.md`](../README.md) for the build modes and the regression
rules, and each case README for its inputs, options and checks.

## Convergence studies

Convergence studies of a case are in `<case>/setup/studies/<study>/`, each with
its own README (the ladder, the options, the quantities checked). The
`summarize.py` next to the studies prints the convergence tables from a folder
of finished run copies:

```bash
python3 tutorials/manufacturedSolutions/bidomain/cube/setup/studies/summarize.py <results folder>
```

The run copies and the tables are generated output and are not tracked. The
automated coverage of a case is its `regression/` directory.

## Bulk, boundary and gradient rows

Every field verifier (monodomain, bidomain, bath bidomain, eikonal) adds these
rows to its error summary for each potential `q` it reports, after its own rows:

| row | quantity |
|---|---|
| `q_bulk`, `q_boundary` | the error of `q` over the cells without, and with, a face on a physical patch of the mesh `q` lives on |
| `gradQ`, `gradQ_bulk`, `gradQ_boundary` | the error of the cell gradient of the solved `q`, with the solver's gradient scheme and boundary values |
| `gradReconQ`, `gradReconQ_bulk`, `gradReconQ_boundary` | the error of the cell gradient of the exact `q`, with exact cell and boundary-face values, through the same scheme: the reconstruction alone |

The physical patches are the sealed wall on a heart sub-mesh, the walls of a
whole-mesh case and the outer patches of a heart-and-bath mesh; processor,
cyclic and empty patches are not. L1 and L2 of every part are normalised by the
total volume of the mesh, so L1 and L2² of the parts add up to the total;
Linf is the maximum over the part. Gauge-shifted potentials (`phiE_gauge`)
keep their field rows under that name and their gradient rows under `phiE`.

On a heart-and-bath mesh the extracellular potential has a third part,
`phiE_interface`, the cells next to a heart–bath face. Its solved gradient is
the one the bath operator consumes (`phiEInterfaceGradient`: jump-aware under
`conormalHarmonic`), and the reconstruction rows leave the interface band out,
since the exact potential has a normal-derivative jump there.

Comparing `gradQ_boundary` with `gradReconQ_boundary` separates an
inconsistent boundary value of the solved field from the accuracy of the
reconstruction next to the boundary.
