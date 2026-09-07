# Distributed 1-D Graph Parallelization Plan

**Status:** proposed  
**Scope:** `conductionSystemDomain` with `monodomain1DSolver`  
**Primary goal:** reduce 1-D graph wall time and per-rank memory while preserving
the existing finite-volume cable equation, reaction/diffusion ordering, PVJ
coupling semantics, dictionary contract, fields, and output artifacts.

## 1. Motivation and verified baseline

The current implementation distributes the ionic ODE work but replicates the
graph diffusion solve.

- `conductionSystemDomain::initialiseState()` assigns each rank a contiguous
  graph-node interval and constructs its ionic model with `nLocalNodes_`
  integration points (`src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.C:260`).
- `monodomain1DSolver::advance()` advances only those local ODE states, places
  the resulting current in a graph-sized array, and reconstructs the full
  `Iion` field with a world-communicator reduction
  (`src/electroModels/conductionSystemModels/monodomain1DSolver/monodomain1DSolver.C:47`).
- Every rank then allocates graph-sized Hines arrays, assembles all edges, and
  performs the complete upward and downward sweeps
  (`src/electroModels/conductionSystemModels/monodomain1DSolver/monodomain1DSolver.C:90`,
  `:178`).
- `Vm1D_`, `Iion1D_`, and `activationTime_` are complete graph fields on every
  rank (`src/electroModels/electroDomains/conductionSystemDomain/conductionSystemDomain.H:80`).
- Tissue-side PVJ voltage gathering uses MPI reductions, while PVJ source
  deposition writes only local tissue cells
  (`src/electroModels/electroCouplers/pvjCoupler/pvjMapper.C:188`, `:243`).

For a tree, `E=N-1`. The replicated diffusion therefore consumes approximately
`O(N)` time and scratch space per rank, or `O(PN)` aggregate diffusion work,
while the current full-field ionic reduction communicates `O(N)` values per
timestep. These are scaling deductions from the loops and storage above, not
measured performance results.

## 2. Scientific and compatibility invariants

Unless separately approved, implementation work must preserve all of the
following:

1. The graph remains a connected tree and continues to reject unsupported
   topology (`conductionGraph.H:146`).
2. The node control length remains the sum of half-lengths of incident edges.
3. Edge flux and the `chi*Cm*controlLength` scaling remain unchanged.
4. The diffusion step remains backward-Euler implicit.
5. Ionic ODE advancement, applied current, PVJ-current sign, and staggered
   graph-before-myocardium ordering remain unchanged.
6. `Vm`, `ionicCurrent`, `activationTime`, terminal currents, and terminal
   sources retain their current units and meanings.
7. Existing dictionary keys and defaults remain valid; new parallel controls
   must have defaults reproducing current behavior.
8. Restart fields, time-series output, VTK output, verification hooks, and PVJ
   coupling remain available.
9. Serial results remain the numerical reference. MPI changes may alter only
   floating-point reduction order within an agreed tolerance.

Any proposed change to an equation, numerical scheme, threshold interpolation,
PVJ kernel normalization, physical constant, unit, or field meaning is outside
this plan and requires a separate scientific decision.

## 3. Design principles

- Establish measurements before selecting a backend.
- Separate graph topology, stimulus placement, ownership, and communication.
- Introduce one behavior boundary at a time and retain a reference backend.
- Partition the tree for the solver first; treat tissue/PVJ locality as a
  secondary partition objective.
- Communicate subtree-interface data rather than complete graph fields.
- Keep all ranks on matched collective paths until communicator splitting is a
  separately validated feature.
- Support the OpenFOAM versions accepted by `Allwmake`; do not assume that a
  communicator API available in one distribution is portable to all of them.

## 4. Phase 0 — correctness blockers and characterization

### 4.1 Separate topology root from stimulus node

Current topology construction always starts at node 0
(`conductionGraph.H:171`), while `rootNode_` is read afterward and can then be
overwritten by `rootStimulus.node`
(`conductionSystemDomain.C:174`, `:248`). The override is not range-checked
after assignment.

Work:

- Add an explicit, validated topology-root input to `conductionGraph` topology
  construction.
- Store topology root and stimulus node as separate domain members.
- Preserve the current dictionary behavior: graph `rootNode` supplies the
  topology root and the default stimulus node; `rootStimulus.node` overrides
  only stimulus placement.
- Validate both values after all dictionary reads and before either is used.
- Add nonzero-root, invalid-root, and independent-stimulus-node tests.

Acceptance:

- Rooting a tree at a different valid node changes only parent/order metadata,
  not the serial solution.
- Invalid graph roots and stimulus nodes fail deterministically with a clear
  diagnostic before state arrays are indexed.

### 4.2 Correct the existing local partition formula

The current integer division assigns the entire graph to the last rank when
`N < P` (`conductionSystemDomain.C:262`). Replace it with a balanced contiguous
partition, for example boundaries `floor(rank*N/P)` and
`floor((rank+1)*N/P)`.

Acceptance:

- Local counts differ by at most one.
- Every global node has exactly one ODE owner.
- `N=0`, `N=1`, `N<P`, non-divisible `N/P`, and serial cases are covered.

### 4.3 Fix parallel ionic export

`conductionSystemDomain::write()` currently loops over all graph nodes while
indexing rank-local ionic state storage (`conductionSystemDomain.C:792`).

Work:

- Fill only the owned global segment from local ionic state.
- Gather or reduce the completed export field before master-only output.
- Apply the same rule to algebraic variables.

Acceptance:

- MPI ionic export has no out-of-bounds access.
- Serial and MPI output contain the same graph-node ordering and values within
  tolerance.

### 4.4 Add characterization tests

Capture before-change trajectories for:

- full graph `Vm` and `Iion`;
- graph activation times;
- terminal voltage, current, and volumetric source;
- selected ionic states;
- restart/read-back behavior;
- serial and existing six-rank coupled Niederer runs.

Comparison must occur at every output time, not only at the final state.

## 5. Phase 1 — measurement and local optimization

### 5.1 Add detailed timing

Add cumulative and per-step timings for:

- local Vm packing;
- ionic ODE solve;
- ionic-current collective;
- applied-current assembly;
- static/dynamic coefficient assembly;
- upward Hines elimination;
- downward substitution;
- activation update;
- PVJ terminal gather and source deposition.

Report minimum, mean, and maximum rank time. The maximum is the timestep
critical path; aggregate CPU time alone must not be used to claim speedup.

### 5.2 Cache graph-invariant data

`controlLength`, edge orientation, and geometry-dependent coefficient factors
are rebuilt every timestep (`monodomain1DSolver.C:90-165`). Cache:

- control length per node;
- child-to-parent edge mapping;
- endpoint coefficient bases excluding `dt` if timesteps may vary;
- reusable scratch buffers.

The diagonal and RHS must still be reset each timestep, and `dt`-dependent
coefficients must be refreshed when `dt` changes.

Acceptance:

- Bitwise-equal serial output where operation ordering is unchanged; otherwise
  agreement within the established tolerance.
- No repeated graph-sized allocations in the timed advance path.
- A measured improvement on at least one representative realistic graph, with
  no statistically significant regression on the small regression graph.

## 6. Phase 2 — selectable master reference backend

Implement a runtime-selectable backend while retaining the replicated solver:

- `replicatedHines`: current behavior and compatibility reference.
- `masterHines`: rank 0 assembles and solves diffusion, followed by full `Vm`
  synchronization.

For `masterHines`:

1. Keep distributed ODE advancement and the current `Iion` collective.
2. Guard all diffusion assembly and both sweeps on master.
3. Broadcast `Vm` before post-graph PVJ coupling.
4. Either update activation on every rank after the `Vm` broadcast or update on
   master and broadcast `activationTime`.
5. Keep every rank entering the same collective calls.

This backend is a correctness oracle and possible resource-saving mode. It is
not considered a scalable speedup unless timing shows lower maximum-rank graph
time than `replicatedHines`; concurrent replicated sweeps may be faster than a
serial solve plus broadcast.

Acceptance:

- `replicatedHines` remains the initial default.
- Serial, replicated MPI, and master MPI trajectories agree.
- Deadlock-free operation for `P=1,2,6` and `N<P`.
- Timings separately expose solve and broadcast costs.

## 7. Phase 3 — explicit graph ownership layer

Introduce an ownership abstraction independently of distributed Hines.

Required data:

- owner rank for each node and edge;
- ordered owned global-node IDs per rank;
- global-to-local and local-to-global maps;
- ghost/interface nodes and their owners;
- terminal-node owners;
- partition adjacency and parent/child relationships;
- send/receive schedules built once after graph construction.

Initially retain global topology metadata if that simplifies validation, but
move evolving state and scratch storage toward owned plus ghost values.

Replace interfaces that assume a single contiguous `localStartNode`, including:

- local Vm packing in `monodomain1DSolver`;
- `manufacturedGraphVerifier`, which currently computes
  `globalNode = localStartNode + localNode`
  (`src/verificationModels/monodomainVerification/manufacturedGraphVerifier.C:143`);
- coupled manufactured-source setup
  (`src/verificationModels/coupledVerification/coupled1D3DMonodomainVerifier.C:362`);
- ionic state and algebraic output.

Communication progression:

1. Preserve replicated `Vm` temporarily while introducing ownership maps.
2. Replace the graph-sized `Iion` all-reduce with ownership-targeted exchange or
   a gather required only by the selected reference backend.
3. Replicate only the small terminal interface where useful.
4. Move `Vm` and activation to owned/ghost storage after the distributed solve
   is validated.

Acceptance:

- Exactly one owner per node and edge.
- Complete, duplicate-free global reconstruction for diagnostics and output.
- No dependence on contiguous global node numbering outside the optional
  baseline partitioner.
- Per-rank evolving graph storage approaches `O(N/P + Nghost + Nterminal)`.

## 8. Phase 4 — distributed Hines solver

### 8.1 Partitioning

Root the validated tree at `topologyRootNode` and divide it into connected
subtrees. The first partitioner should be deterministic and internal to the
repository; an external graph-partitioning dependency is not required for the
initial implementation.

Objectives, in order:

1. balance owned nodes and estimated ionic-model cost;
2. minimize cut edges;
3. preserve useful branch-level concurrency;
4. optionally improve PVJ/tissue locality.

The partition quotient must itself be a tree. If one rank owns multiple
subtrees, each component must be represented explicitly rather than hidden by
a single-range assumption.

### 8.2 Upward elimination

For each owned subtree:

- assemble local diagonal, RHS, and parent/child coefficients;
- eliminate strictly interior descendants locally;
- condense each cut edge to the interface relation required by its parent;
- send the condensed relation toward the parent-owner rank;
- continue until the partition containing the topology root can solve the root
  interface value.

Batch all relations sent to the same rank and communication level. Avoid a
blocking message for every graph node.

### 8.3 Downward substitution

- Send solved parent-interface voltages to child-subtree owners.
- Perform local top-down substitution once the required parent value arrives.
- Exchange only the ghost/interface values required by dependent subtrees and
  terminal consumers.

### 8.4 Numerical safeguards

- Preserve the existing asymmetric endpoint coefficients.
- Detect zero or invalid pivots before division and report the global node and
  owner rank.
- Use deterministic traversal and packing order.
- Keep the master backend available for step-by-step comparison.
- Add optional debug reconstruction of full matrix residual
  `||A*Vm-rhs||` for small graphs.

Acceptance:

- Distributed and master residuals meet the same tolerance.
- Full trajectories agree for balanced, unbalanced, deep-chain, star, and
  realistic branching trees.
- Communication volume scales with partition interfaces rather than all graph
  nodes, excluding requested global output times.
- Measured wall-clock speedup is reported against both replicated and master
  backends.

## 9. Phase 5 — PVJ communication optimization

The current mapper performs one reduction per terminal
(`pvjMapper.C:197-207`). Once graph ownership is stable:

- pack all terminal tissue contributions into one buffer and reduce it in one
  collective;
- route network terminal voltage only to ranks that need it;
- retain local tissue deposition;
- make terminal ownership and any replicated terminal cache explicit;
- verify conservation of current deposited into all local tissue cells.

Graph-to-tissue affinity may then be added as a secondary partition weight.
Do not force mesh-coincident ownership if it materially increases tree cut
edges or the distributed-Hines critical path.

Acceptance:

- Terminal gather uses a bounded number of collectives per coupling phase,
  independent of terminal count.
- Terminal voltage/current and total deposited source match the reference.

## 10. Phase 6 — optional dedicated graph communicator

Consider graph-only ranks only if profiling demonstrates that graph work is a
material fraction of total runtime on realistic cases and phases 1-5 are
insufficient.

Required before implementation:

- a portable communicator wrapper validated against every supported OpenFOAM
  version family;
- explicit world-to-graph routing for tissue terminal voltage;
- graph-to-world routing for terminal voltage/current/source and activation;
- collective failure propagation;
- construction, restart, write, and shutdown behavior for ranks outside the
  graph communicator;
- load-balance analysis showing that reserving graph ranks does not slow the 3-D
  solve more than it accelerates the graph.

Current orchestration calls coupling and graph phases on all ranks
(`staggeredElectrophysicsAdvanceScheme.C:66-75`), and the current ionic and PVJ
collectives are world-wide. Consequently, communicator splitting is an
architectural phase, not a rank guard around existing code.

## 11. Test and benchmark matrix

### 11.1 Correctness matrix

Run each applicable backend with `P=1,2,6` and at least one `P>N` case:

| Case | Purpose | Required comparisons |
| --- | --- | --- |
| Two-node tree | Minimal edge and root behavior | full fields, residual |
| Linear chain | Worst-case tree depth | full trajectory, deadlock |
| Star tree | Maximum branch concurrency | interfaces, residual |
| Uneven branching tree | Partition balance | ownership, full trajectory |
| Nonzero topology root | Root correctness | parent/order, full trajectory |
| Separate stimulus node | Root/stimulus semantics | applied current, Vm |
| Manufactured 1-D/1-D–3-D | Accuracy | error norms and convergence order |
| Niederer Purkinje slab | Coupled regression | graph/PVJ reference values |
| Restart case | Artifact compatibility | continuous vs restarted trajectory |
| Ionic export enabled | MPI output correctness | all state/algebraic columns |

The manufactured 1-D–3-D suite documents expected second-order spatial
convergence (`tutorials/manufacturedSolutions/monodomain1D3D/README.md:87`).
The Niederer regression already performs a six-rank coupled run and checks
graph/PVJ values
(`tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/regression/regressionTest.sh:223`).

### 11.2 Performance matrix

Measure at minimum:

- graph sizes spanning small regression, medium, and realistic PATHOS-scale
  networks;
- ranks `P=1,2,4,6,8,...` up to the useful graph concurrency;
- replicated, master, and distributed backends;
- graph-only timed kernels and complete coupled timestep time;
- maximum-rank time, aggregate CPU time, communication time, and peak resident
  memory per rank;
- strong-scaling efficiency and time spent idle or waiting.

Use multiple repetitions after warm-up. Record compiler, OpenFOAM version,
hardware, MPI implementation, graph characteristics, timestep, ionic model,
and output settings.

### 11.3 Performance gates

- No backend becomes the default solely because it reduces aggregate work.
- `masterHines` becomes preferred only where it reduces measured critical-path
  time or solves a documented memory/resource constraint.
- `distributedHines` becomes the default only after correctness passes on all
  supported cases and it shows a repeatable coupled-runtime benefit on target
  workloads.
- If graph time is negligible relative to the myocardium solve, retain the
  simplest correct backend and avoid communicator complexity.

## 12. Proposed code organization

Keep solver policy out of `conductionSystemDomain` where possible:

- `conductionGraph`: validated topology and root-oriented traversal.
- `graphPartition`: ownership, local/global maps, ghosts, and schedules.
- `graphLinearSolver` interface: replicated, master, and distributed Hines
  backends.
- `monodomain1DSolver`: reaction/diffusion orchestration and physical assembly.
- `conductionSystemDomain`: state, dictionary contract, coupling endpoints,
  restart, and output.
- `pvjMapper`/PVJ couplers: tissue mapping and terminal exchange.

Introduce backend selection only after the reference implementation exists.
The default must initially remain current replicated behavior.

## 13. Delivery batches

Each batch must be independently reviewable and revertible:

1. Root/stimulus separation and validation tests.
2. Balanced contiguous ODE partition and `N<P` tests.
3. Parallel ionic-export correction.
4. Fine-grained timing and baseline report.
5. Static coefficient caching and scratch reuse.
6. Selectable master backend and MPI equivalence tests.
7. Ownership/maps abstraction with current numerical backend.
8. Ownership-aware ionic-current communication.
9. Distributed Hines upward sweep.
10. Distributed Hines downward sweep and state localization.
11. Packed PVJ terminal communication.
12. Default-backend decision based on benchmark evidence.
13. Optional communicator prototype, gated by profiling.

Do not combine a scientific-model change with any of these batches.

## 14. Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| Deadlock from mismatched collectives | Centralize schedules; test empty ownership and `N<P`; add phase-tagged diagnostics |
| Loss of numerical reproducibility | Deterministic traversal/packing; compare every output time with master backend |
| Root/stimulus semantic regression | Separate members and dedicated nonzero-root tests |
| Output/restart incompatibility | Reconstruct global ordering only at artifact boundaries; retain existing field names |
| Poor scaling on deep trees | Benchmark chain topology; consider interface-level parallel algorithms only if required |
| Too many small MPI messages | Batch by neighbor and traversal level; precompute schedules |
| PVJ locality degrades Hines partition | Make tissue affinity a weighted secondary objective |
| OpenFOAM communicator API divergence | Isolate MPI/UPstream calls behind a small compatibility layer and compile every supported version |
| Optimization is irrelevant to coupled runtime | Require end-to-end timings and a performance gate before changing defaults |

## 15. Definition of done

The project is complete when:

1. root and stimulus semantics are explicit and validated;
2. serial, replicated, master, and distributed results pass the agreed
   equivalence and manufactured-accuracy tests;
3. MPI ionic output and restart behavior are correct;
4. distributed graph state and scratch memory scale with owned plus ghost data;
5. steady-state timestep communication is proportional to partition interfaces
   and terminal exchange rather than full graph fields;
6. realistic coupled benchmarks demonstrate a repeatable wall-clock benefit;
7. the selected default backend is justified by recorded evidence;
8. all existing dictionary keys, units, fields, coupling behavior, and output
   artifacts remain compatible.

