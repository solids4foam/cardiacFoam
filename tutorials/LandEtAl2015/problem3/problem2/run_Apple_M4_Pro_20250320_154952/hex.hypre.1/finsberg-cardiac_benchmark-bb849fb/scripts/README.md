# Scripts

This folder contains script for submitting the jobs on the [EX3 cluster](https://www.ex3.simula.no). However, the scripts should be fairly generic since they use SLURM. The only thing you probably need to change in order to run on a different clusters are the modules that you load in the beginning and the file names.

## Submit benchmark

Step 0 (case A + B)
```
python run_benchmark.py step0
```

Step 1
```
python run_benchmark.py step1
```

Step 2 (case A + B + C)
```
python run_benchmark.py step2
```

Run all benchmarks
```
python run_benchmark.py all
```

Note that it is also possible to add a the flag `--dry-run` to just check which benchmarks that are run.
