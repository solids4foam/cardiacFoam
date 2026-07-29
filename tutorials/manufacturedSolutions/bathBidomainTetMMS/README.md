# manufacturedSolutions/bathBidomainTetMMS

For the selected formulation, final four-level results, limitations, and
paper-ready conclusion, see [`FINAL_SOLUTION.md`](FINAL_SOLUTION.md). The files
under `setup/interfaceStudy/` retain the underlying verification audit trail.

Tetrahedral-mesh extension of the structured `bathBidomain` manufactured
solution. The mesh gate is followed by a separate N=10 solver smoke test.

The geometry matches the structured reference over `[-1,2] x [0,1] x [0,1]`:

- left bath: `-1 <= x <= 0`
- myocardium: `0 <= x <= 1`
- right bath: `1 <= x <= 2`

All three volumes are constructed and fragmented in one Gmsh model. The two
heart--bath interfaces remain internal conformal faces. The imported OpenFOAM
mesh has two cellZones, `myocardium` and `bath`, and three exterior patches,
`xMin`, `xMax`, and `sides`.

Run the N=10 mesh gate with:

```bash
bash setup/run_mesh_gate.sh 10
```

Results are archived under `setup/results/N10/`. Review `mesh_manifest.txt`,
`log.checkMesh`, `log.checkMesh.strict`, and `log.gmshToFoam` before adding
solver dictionaries or starting a convergence sweep. The standard mesh check
is the pass/fail gate; the extended topology/geometry audit is also retained
to expose low-determinant tetrahedral boundary stencils.

Run the complete N=10 solver smoke test with:

```bash
bash Allrun.smoke
```

The smoke run uses the structured bath-bidomain manufactured problem unchanged
apart from `dimension "3D"`, tetrahedral spatial schemes, and the shortened
N=10 integration window. Solver logs and manufactured summaries are archived
under `setup/results/N10/smoke/`.

After the smoke test passes, run the three-point spatial sweep with:

```bash
bash setup/run_tet_sweep.sh
```

The default ladder is `N=10 20 40`, with the structured-case `deltaT ~ h^2`
schedule and `endTime=0.02`. Each resolution is generated from scratch and
archived under `setup/results/N<N>/`. The combined `summary.csv` computes
observed order from `h_heart=(1/N_myocardium_cells)^(1/3)`, not from nominal
Gmsh `N` or the verifier's rounded structured-grid filename.

`run_tet_sweep.sh` reproduces the **potential** convergence half of the paper's
tetrahedral bath-bidomain table (`Vm`, `phiE`, `phiI` errors -> `summary.csv`).
The **assembled-current / interface-flux** half (the interface-conservation
result) comes from the `bathBidomainInterfaceMetrics` utility, which
`run_tet_sweep.sh` does not run. Reproduce that half with:

```bash
ASSEMBLY=matchedSubmesh METHODS=distanceWeightedHarmonic RESOLUTIONS="10 20 40 80" \
  bash setup/run_parallel_interface_sweep.sh
```

which writes `setup/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic/N<N>/bathBidomainInterfaceMetrics.csv`
(the source of record; see `FINAL_SOLUTION.md`). Both halves use the same
selected formulation (`matchedSubmesh` + `distanceWeightedHarmonic`, baked into
`constant/electroProperties`) and the committed `snGrad corrected` scheme.
