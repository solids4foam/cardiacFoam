# GPU performance characterization

## Scope and environment

This is a measured characterization of CUDA active tension, not a claim that
the entire monodomain PDE or mechanics solver is GPU-resident. The workload is
the 3,360-cell `springSupportedSlab` with a 10 us shared coupled step through
20 ms (2,000 updates), one MPI rank, OpenFOAM v2412, and the isolated target
libraries in `/tmp/cardiac-gpu-characterization-20261001/lib`. The node GPU was
an NVIDIA RTX 4000 Ada Generation (20,475 MiB), driver 575.57.08.

The case's TNNP ionic model, tissue solve, and solids4foam mechanics remain
host-resident. Thus wall time measures end-to-end impact in the real coupled
case, while CUDA events isolate active-tension work.

## End-to-end A/B timing

| Active-tension model | Host-batched [s] | CUDA [s] | CUDA speedup | Result |
| --- | ---: | ---: | ---: | --- |
| LandNiedererTWorldBatched | 149.43 | 139.64 | 1.070x | identical 2,000-sample `Ta` probe trace |
| NashPanfilovBatched | 153.98 | 144.38 | 1.066x | identical 2,000-sample `Ta` probe trace |

These are unprofiled complete solver wall times. They include initialization,
host PDE and mechanics work, I/O, and active-tension execution. They are not
pure ODE-kernel throughput measurements.

## TWorld CUDA event profile

`batchedCUDAProfile true` was enabled only for a fresh TWorld CUDA run. The
instrumentation records steady-state model inputs, all ODE kernels (including
the final post-step RHS evaluation), and the returned `Ta` row. It excludes
one-time state allocation/upload from H2D, but includes it in wrapper total.

| Component | Aggregate [s] | Per 10-us update [us] | Share of 141.00 s profiled case |
| --- | ---: | ---: | ---: |
| H2D inputs | 0.103586 | 51.8 | 0.073% |
| ODE kernels | 0.147324 | 73.7 | 0.104% |
| D2H `Ta` | 0.0368411 | 18.4 | 0.026% |
| Named device operations | 0.287751 | 143.9 | 0.204% |
| Active-tension wrapper | 0.542988 | 271.5 | 0.385% |

The profile uses CUDA events and synchronizes each timed segment, so it
perturbs execution. Use unprofiled A/B rows for whole-case timing and event
values only to rank bottlenecks. At this mesh size CPU tissue/mechanics
dominate; immediate GPU targets are reducing the three input rows per update
and kernel-launch aggregation, after maintaining trace parity.

## Reproduction

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export LD_LIBRARY_PATH=/tmp/cardiac-gpu-characterization-20261001/lib:$LD_LIBRARY_PATH

# In a copy of tutorials/electromechanicsProtocols/springSupportedSlab:
# select LandNiedererTWorldBatched and set batchedCUDAProfile true.
blockMesh
cp -r constant/polyMesh constant/electro/polyMesh
cp -r constant/polyMesh constant/solid/polyMesh
setExprFields -region electro -time 0 -ascii
srun -N1 -n1 -p dev cardiacFoam | tee log.cardiacFoam
```

The shutdown line has the form `CUDA profile: calls=... H2D_s=...`.
Disable `batchedCUDAProfile` for normal performance measurements.

## Limits and next measurements

- This does not characterize existing ionic CUDA kernels. Their device
  residency and transfers must be measured separately on myocardial slabs and
  Niederer geometry.
- The direct GPU benefit is modest because only active-tension ODEs are
  offloaded in this case; it is not full-solver GPU speedup.
- Long deforming LandNiederer-2017 spring runs are currently invalid because
  scalar and host-batched implementations both fail nonlinear mechanics
  convergence near 18 ms. Do not use this model in end-to-end performance
  comparisons until the coupling issue is resolved.
