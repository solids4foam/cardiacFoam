# MeluXina Benchmark Plan

This folder contains MeluXina-oriented benchmark planning and Slurm templates
for the slab cases under `tutorials/benchmarkGPU/cases`.

Reference docs used:

- System overview: GPU accelerator nodes provide `4x NVIDIA Ampere 40GB HBM`
  per node, CPU cluster nodes provide high core-count CPU-only nodes.
- Job handling: use Slurm `--account`, `--partition`, `--qos`, and for GPU/MPI
  jobs start from `--gpus-per-task=1` and one MPI rank per GPU.

Links:

- <https://docs.lxp.lu/system/overview/>
- <https://docs.lxp.lu/first-steps/handling_jobs/>

## Resource model

- CPU partition: `cpu`
- GPU partition: `gpu`
- Default QOS for production-style runs: `default`
- Test QOS for short smoke tests: `test`

MeluXina hardware relevant here:

- CPU nodes: 128 physical CPU cores per node
- GPU nodes: 4 GPUs per node

That implies these starting mappings:

- CPU `np=1,8,24,48` all fit on one CPU node
- GPU/MPI `np=8` => 2 GPU nodes
- GPU/MPI `np=24` => 6 GPU nodes
- GPU/MPI `np=48` => 12 GPU nodes

## Recommended run matrix

### CPU adaptive reference runs

| Model | Case | MPI ranks | Node type | Nodes |
|---|---|---:|---|---:|
| BuenoOrovio | `cpu_bueno_adaptive_dt2e-5` | 1, 8, 24, 48 | CPU | 1 |
| TNNP | `cpu_tnnp_adaptive_dt2e-6` | 1, 8, 24, 48 | CPU | 1 |
| TWorld | `cpu_tworld_adaptive_dt2e-6` | 1, 8, 24, 48 | CPU | 1 |

### GPU serial runs

| Model | Case | MPI ranks | GPU nodes | GPUs total |
|---|---|---:|---:|---:|
| BuenoOrovio | `gpu_bueno_compact_rl2_dt2e-5` | 1 | 1 | 1 |
| TNNP | `gpu_tnnp_compact_rl10_dt2e-5` | 1 | 1 | 1 |
| TWorld | `gpu_tworld_compact_rl10_dt2e-5` | 1 | 1 | 1 |

### GPU+MPI runs

These are the compact hybrid points that looked most informative locally.

| Model | Case | MPI ranks | GPU nodes | GPUs total |
|---|---|---:|---:|---:|
| BuenoOrovio | `gpu_bueno_compact_rl2_dt2e-5` | 8, 24, 48 | 2, 6, 12 | 8, 24, 48 |
| TNNP | `gpu_tnnp_compact_rl10_dt1e-5` | 8, 24, 48 | 2, 6, 12 | 8, 24, 48 |
| TWorld | `gpu_tworld_compact_rl10_dt1e-5` | 8, 24, 48 | 2, 6, 12 | 8, 24, 48 |

## What you need to prepare on MeluXina

Before submitting:

1. Copy the repo or benchmark case directory to your MeluXina project space.
2. Build `cardiacFoam` in the MeluXina software environment you intend to use.
3. Identify:
   - your Slurm account, e.g. `p20xxxx`
   - target QOS: `default` or `test`
4. Check module environment for:
   - compiler / MPI
   - CUDA or NVIDIA HPC stack if needed
5. Verify the binary and GPU library path on a short interactive test.

Minimum command checks on MeluXina:

```bash
sacctmgr show user $USER withassoc format=user,account,defaultaccount
scontrol show partition cpu
scontrol show partition gpu
```

## Templates

- `meluxina_cpu_case.slurm`
- `meluxina_gpu_case.slurm`
- `meluxina_gpu_mpi_case.slurm`

All templates assume you submit them from `tutorials/benchmarkGPU/meluxina/`
and pass the benchmark case name through `BENCHMARKGPU_CASE`.

## Suggested first executions

If you want the smallest sensible MeluXina validation set first:

1. CPU:
   - `cpu_bueno_adaptive_dt2e-5` at `np=1`
   - `cpu_tnnp_adaptive_dt2e-6` at `np=1`
   - `cpu_tworld_adaptive_dt2e-6` at `np=1`
2. GPU serial:
   - `gpu_bueno_compact_rl2_dt2e-5`
   - `gpu_tnnp_compact_rl10_dt2e-5`
   - `gpu_tworld_compact_rl10_dt2e-5`
3. GPU+MPI:
   - `gpu_bueno_compact_rl2_dt2e-5` at `np=8`
   - `gpu_tnnp_compact_rl10_dt1e-5` at `np=8`
   - `gpu_tworld_compact_rl10_dt1e-5` at `np=8`

That gives you:

- one CPU reference point per model
- one GPU serial point per model
- one hybrid point per model

before committing to the full `8/24/48` sweep.
