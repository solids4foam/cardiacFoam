# Slab post-processing

These scripts analyze outputs produced by the `ioniGPUKernel` slab runners.
The runners and case generators are kept in the parent tutorial directory;
`setup` contains the post-processing tools used by those runners.

## Field and activation parity

[`post_processing_slab_field_parity.py`](post_processing_slab_field_parity.py)
compares two OpenFOAM case output directories at matching written times. Give
the reference case first and the candidate case second:

```bash
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_slab_field_parity.py \
    /path/to/reference-case /path/to/candidate-case \
    --times 0.005 0.01 0.015
```

For each time, the report contains:

- `Vm_RMSE_mV` and `Vm_max_mV`: candidate minus reference over all cells.
- `activated_scalar` and `activated_batched`: activated counts, where
  `activationTime > 0`.
- `activation_mask_mismatch`: number of cells activated in only one case.
- `activation_p95_ms`: 95th percentile of absolute activation-time differences
  for cells activated in both cases.

The output keys `scalar` and `batched` are historical names required by the
sweep parser. They denote the first and second cases, respectively; they do
not identify the execution backend. The activation-mask mismatch is reported
separately because equal total counts alone do not prove that the same cells
activated. Voltage RMSE and maximum error are unaligned pointwise comparisons,
so a traveling wavefront can produce large values from a small timing shift.

This comparison is appropriate for scalar-versus-batched and
host-batched-versus-CUDA-batched runs when mesh, stimulus, output times, and
case configuration match. It does not itself test ionic-current parity, ODE
state parity, APD, or propagation speed.

## Integrator sweep summary

[`post_processing_integrator_sweep_summary.py`](post_processing_integrator_sweep_summary.py)
reads `summary.csv` from `run_slab_integrator_study.py`, selects the fastest
completed CUDA settings that pass the provisional activation screen, and
writes `best_settings.md` in that sweep output directory:

```bash
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_integrator_sweep_summary.py \
    /path/to/integrator-sweep-output
```

A setting passes only if each of the three configured output times is
available, activated-cell mask mismatch is at most
`max(5 cells, 1% of slab cells)`, and activation-time p95 is at most 0.1 ms
against both the matched-step scalar RKF45 case and the fine-step scalar
reference. The screen does not threshold pointwise voltage error because a
traveling upstroke can dominate that metric. It is a short propagation screen,
not full APD/state/current validation or evidence of production stability.

Reported runtime is compared with the matched-step scalar run and includes the
case `Allrun`, including mesh generation and the CPU tissue solve. Use the
per-case `summary.csv` for failures and detailed comparisons; a selected setting
still needs longer physiological validation before production use.
