# Active-tension CUDA validation

## Scope

These results cover the three active-tension families that have scalar,
batched-host, and CUDA implementations: `NashPanfilov`, `LandNiederer`, and
`LandNiedererTWorld`.  The tissue PDE and mechanics solver remain host
resident; this is GPU acceleration of the batched active-tension ODEs only.

## Test environment

- Source worktree: `gpu-characterization-20261001` from `a6813e26`.
- OpenFOAM: v2412; CUDA compiler: nvcc with GCC 10 host compiler.
- GPU: NVIDIA RTX 4000 Ada, one rank on the local `dev` Slurm partition.
- Active-tension library:
  `/tmp/cardiac-gpu-characterization-20261001/lib/libactiveTensionModels.so`.

## Controlled single-cell parity

`tutorials/electrophysiologyProtocols/singleCell` was run for 100 ms with a
1 us solver step and the same TWorld ionic drive, stimulus, initial state, and
parameters in each variant.  Scalar models use their existing ODE solver;
batched host and CUDA models use the documented explicit Euler update with one
substep.  The exported `Ta` trace contains 101 samples.

| Model | Scalar vs host-batched max \|Ta\| [kPa] | Host-batched vs CUDA | CUDA device dispatch | Result |
| --- | ---: | ---: | --- | --- |
| NashPanfilov | 0.0067594 | 0 at 101/101 samples | yes | pass |
| LandNiederer | 8.04e-05 | 0 at 101/101 samples | yes | pass |
| LandNiedererTWorld | 5.21e-05 | 0 at 100 written samples | yes | pass |

The scalar/batched values are integration-method differences, not CUDA
differences: every batched host/CUDA pair is exactly equal in the written
trace.  The scalar comparison is acceptable for these traces; it must still
be repeated with timestep-convergence and event metrics before being used as a
physiological accuracy claim.

## Coupled spring-supported slab gate

The `springSupportedSlab` case has 3,360 integration cells and a 10 us coupled
step.  The following are short 0.2 ms (20-step) checks, not the documented
0.25 s mechanics benchmark.

| Model | Host-batched | CUDA | Mechanics | Host/CUDA Ta probes | Result |
| --- | --- | --- | --- | --- | --- |
| NashPanfilovBatched | End | End; device logged | converged every step | identical | pass |
| LandNiedererBatched | End | End; device logged | converged every step | identical | pass |
| LandNiedererTWorldBatched | End | End; device logged | converged every step | identical | pass |

This also verifies the corrected Land-family preconditioning path.  The former
100 x 10 ms Euler preconditioner was unstable at TNNP resting calcium; batched
Land models now default to a maximum 0.1 ms startup increment through
`batchedPreconditioningMaxStep`.

## CUDA residency and performance interpretation

CUDA ODE state, rates, and algebraics now remain resident after the first
update.  Per solve, the implementation transfers only current model inputs and
the `Ta` row required by the host coupled solver.  Full state/rate/algebraic
copies are performed only when the existing restart/export I/O interface asks
for them.  A post-change 20-step LandNiedererTWorld CUDA slab reproduced the
host probe file exactly and completed in 2.83 s wall time (including Slurm
launch and complete solver work).

The recorded 100 ms single-cell CUDA wall times—4.86 s (Nash), 6.91 s (Land),
and 6.88 s (LandTWorld)—include solver startup, host ionic work, I/O, and Slurm
launch.  They are not an ionic- or active-tension-kernel speedup measurement.
Kernel-level timings and an A/B transfer benchmark on a representative slab
remain required before reporting acceleration.

## Reproduction

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export LD_LIBRARY_PATH=/tmp/cardiac-gpu-characterization-20261001/lib:$LD_LIBRARY_PATH
blockMesh -case tutorials/electrophysiologyProtocols/singleCell
srun -N1 -n1 -p dev cardiacFoam -case tutorials/electrophysiologyProtocols/singleCell
```

Select `NashPanfilovBatched`, `LandNiedererBatched`, or
`LandNiedererTWorldBatched` in the case dictionary.  Leave `batchedUseCUDA`
at its default for CUDA, or set it to `false` for the host-batched reference.
For the coupled gate, follow `springSupportedSlab/Allrun` after selecting the
same model in `constant/electroMechanicalProperties`.

## Remaining work

- Run the complete 0.25 s coupled benchmark and timestep-convergence study.
- Add GPU-event timing around H2D, ODE kernels, D2H, and synchronization.
- Validate restart continuation across a write/restart boundary on CUDA.
- Complete the requested per-ionic-model single-cell/slab matrix and the
  Niederer benchmark characterization.
