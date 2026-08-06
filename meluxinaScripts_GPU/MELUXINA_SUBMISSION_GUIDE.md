# MeluXina Submission Guide

This guide explains how to run the `tutorials/benchmarkGPU` slab benchmarks on MeluXina.

Related files in this folder:

- `README.md`
- `meluxina_cpu_case.slurm`
- `meluxina_gpu_case.slurm`
- `meluxina_gpu_mpi_case.slurm`

Reference documentation:

- <https://docs.lxp.lu/system/overview/>
- <https://docs.lxp.lu/first-steps/handling_jobs/>

## What Each Run Type Means

### CPU only

- Uses the CPU partition.
- Uses `np = number of MPI ranks`.
- No GPU is requested.

### GPU serial

- Uses the GPU partition.
- Uses `np=1`.
- Uses `1` MPI rank and `1` GPU.

### GPU + MPI

- Uses the GPU partition.
- Uses `1 MPI rank : 1 GPU`.
- So:
  - `np8` means `8` MPI ranks and `8` GPUs
  - `np24` means `24` MPI ranks and `24` GPUs
  - `np48` means `48` MPI ranks and `48` GPUs

On MeluXina GPU nodes, there are `4` GPUs per node, so:

- `np8` => `2` GPU nodes
- `np24` => `6` GPU nodes
- `np48` => `12` GPU nodes

This is different from the local workstation results. Locally, many MPI ranks shared one GPU. On MeluXina, the hybrid runs are meant to scale across multiple GPUs.

## Files To Edit Before Submission

You should edit these templates first:

- [meluxina_cpu_case.slurm](/home/simao/cardiacFoam/tutorials/benchmarkGPU/meluxina/meluxina_cpu_case.slurm:1)
- [meluxina_gpu_case.slurm](/home/simao/cardiacFoam/tutorials/benchmarkGPU/meluxina/meluxina_gpu_case.slurm:1)
- [meluxina_gpu_mpi_case.slurm](/home/simao/cardiacFoam/tutorials/benchmarkGPU/meluxina/meluxina_gpu_mpi_case.slurm:1)

Update at least:

1. `#SBATCH --account=p20xxxx`
2. module loading section
3. walltime if needed
4. optional `--qos` if you want `test` instead of `default` for smoke tests

## Required Environment Checks

Run these first on MeluXina:

```bash
sacctmgr show user $USER withassoc format=user,account,defaultaccount
scontrol show partition cpu
scontrol show partition gpu
```

Also verify:

```bash
which cardiacFoam
ldd $(which cardiacFoam)
nvidia-smi
```

For the GPU jobs, verify that your MeluXina build can load the CUDA-linked ionic model library used by the GPU slab cases.

## Repository Layout Expected By The Templates

The templates assume you submit from:

`tutorials/benchmarkGPU/meluxina/`

and that benchmark cases are under:

`tutorials/benchmarkGPU/cases/<case>`

The job copies the selected case into a private temporary directory, runs there, and copies results back into:

- CPU:
  - `cases/<case>/meluxinaResults/cpu/np<N>/job<JOBID>/`
- GPU serial:
  - `cases/<case>/meluxinaResults/gpu-serial/job<JOBID>/`
- GPU+MPI:
  - `cases/<case>/meluxinaResults/gpu-mpi/np<N>/job<JOBID>/`

## Recommended Cases

### CPU adaptive references

- `cpu_bueno_adaptive_dt2e-5`
- `cpu_tnnp_adaptive_dt2e-6`
- `cpu_tworld_adaptive_dt2e-6`

### GPU serial references

- `gpu_bueno_compact_rl2_dt2e-5`
- `gpu_tnnp_compact_rl10_dt2e-5`
- `gpu_tworld_compact_rl10_dt2e-5`

### GPU + MPI references

- `gpu_bueno_compact_rl2_dt2e-5`
- `gpu_tnnp_compact_rl10_dt1e-5`
- `gpu_tworld_compact_rl10_dt1e-5`

MPI ranks to test:

- `8`
- `24`
- `48`

## First Jobs To Run

Start with a minimal validation set.

### CPU serial checks

```bash
cd /path/to/cardiacFoam/tutorials/benchmarkGPU/meluxina

BENCHMARKGPU_CASE=cpu_bueno_adaptive_dt2e-5 sbatch meluxina_cpu_case.slurm
BENCHMARKGPU_CASE=cpu_tnnp_adaptive_dt2e-6 sbatch meluxina_cpu_case.slurm
BENCHMARKGPU_CASE=cpu_tworld_adaptive_dt2e-6 sbatch meluxina_cpu_case.slurm
```

### GPU serial checks

```bash
cd /path/to/cardiacFoam/tutorials/benchmarkGPU/meluxina

BENCHMARKGPU_CASE=gpu_bueno_compact_rl2_dt2e-5 sbatch meluxina_gpu_case.slurm
BENCHMARKGPU_CASE=gpu_tnnp_compact_rl10_dt2e-5 sbatch meluxina_gpu_case.slurm
BENCHMARKGPU_CASE=gpu_tworld_compact_rl10_dt2e-5 sbatch meluxina_gpu_case.slurm
```

### GPU + MPI smoke tests

```bash
cd /path/to/cardiacFoam/tutorials/benchmarkGPU/meluxina

BENCHMARKGPU_CASE=gpu_bueno_compact_rl2_dt2e-5 \
sbatch --nodes=2 --ntasks=8 --ntasks-per-node=4 meluxina_gpu_mpi_case.slurm

BENCHMARKGPU_CASE=gpu_tnnp_compact_rl10_dt1e-5 \
sbatch --nodes=2 --ntasks=8 --ntasks-per-node=4 meluxina_gpu_mpi_case.slurm

BENCHMARKGPU_CASE=gpu_tworld_compact_rl10_dt1e-5 \
sbatch --nodes=2 --ntasks=8 --ntasks-per-node=4 meluxina_gpu_mpi_case.slurm
```

## Full Hybrid Submit Examples

### `np24`

```bash
BENCHMARKGPU_CASE=gpu_tnnp_compact_rl10_dt1e-5 \
sbatch --nodes=6 --ntasks=24 --ntasks-per-node=4 meluxina_gpu_mpi_case.slurm
```

### `np48`

```bash
BENCHMARKGPU_CASE=gpu_tnnp_compact_rl10_dt1e-5 \
sbatch --nodes=12 --ntasks=48 --ntasks-per-node=4 meluxina_gpu_mpi_case.slurm
```

## Important MPI Launcher Note

MeluXina documentation recommends `srun` for MPI job steps.

Your benchmark `Allrun` files currently rely on OpenFOAM `runParallel`, for example:

- [cpu_bueno_adaptive_dt2e-5/Allrun](/home/simao/cardiacFoam/tutorials/benchmarkGPU/cases/cpu_bueno_adaptive_dt2e-5/Allrun:1)

That means the first real check on MeluXina should be:

1. submit one short MPI job
2. inspect the log
3. confirm that the parallel launch is clean in that environment

If `runParallel` resolves properly under MeluXina's OpenFOAM environment, the templates can be used as-is. If not, the only likely patch needed is to adapt the parallel launch path to the site-preferred `srun`.

## Suggested Execution Order

1. CPU `np1` for the three models
2. GPU serial for the three models
3. GPU+MPI `np8` for the three models
4. If those are clean, extend to `np24`
5. Then extend to `np48`

That gives you:

- basic CPU reference validation
- basic GPU validation
- first multi-GPU scaling point
- then the larger hybrid scaling sweep

## Local Analysis Files To Compare Against

Use these local benchmark summaries when comparing MeluXina results:

- [BENCHMARK_ANALYSIS_2026-06-28.md](/home/simao/cardiacFoam/tutorials/benchmarkGPU/BENCHMARK_ANALYSIS_2026-06-28.md:1)
- [benchmarkGPU_results_summary.csv](/home/simao/cardiacFoam/tutorials/benchmarkGPU/benchmarkGPU_results_summary.csv:1)
- [benchmarkGPU_results_summary.json](/home/simao/cardiacFoam/tutorials/benchmarkGPU/benchmarkGPU_results_summary.json:1)
