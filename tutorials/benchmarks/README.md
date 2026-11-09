# Benchmarks

Cases that reproduce a published benchmark, with the quantities the benchmark
reports.

| Case | Benchmark |
| --- | --- |
| [`NiedererEtAl2011verification`](NiedererEtAl2011verification/README.md) | the Niederer et al. (2011) monodomain slab: activation times and smoke-check fields |

```bash
cd tutorials/benchmarks/NiedererEtAl2011verification
./Allrun
```

The regression, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight \
    ./tutorials/Alltest-regression --only 'benchmarks/NiedererEtAl2011verification'
```
