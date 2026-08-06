#!/bin/bash
driverFoam sweep-run --spec setup/studies/temporalConvergence/sweep_temporal_convergence.json --output-dir setup/studies/temporalConvergence/results/sweepCases --max-cases 1000
