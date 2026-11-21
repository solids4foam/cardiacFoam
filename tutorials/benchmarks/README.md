# Benchmarks

Cases that reproduce a published benchmark, with the quantities the benchmark
reports.

| Case | Benchmark |
| --- | --- |
| [`NiedererEtAl2011verification`](NiedererEtAl2011verification/README.md) | the Niederer et al. (2011) monodomain slab: activation times and smoke-check fields |
| [`NiedererEtAl2011verification/bidomain`](NiedererEtAl2011verification/bidomain/README.md) | the same slab with the bidomain model, with the published conductivities or an equal-anisotropy set that reproduces the monodomain case |

```bash
cd tutorials/benchmarks/NiedererEtAl2011verification
./Allrun
```

The regression is in the `full` tier; from the repository root:

```bash
CARDIAC_REGRESSION_SCOPE=full CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'benchmarks/NiedererEtAl2011verification'
```
