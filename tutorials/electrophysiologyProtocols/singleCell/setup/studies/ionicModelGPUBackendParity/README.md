# Ionic-model GPU backend parity — single cell

This study is the 0-D companion to
[`ionicModelGPUBackendParity`](../../../ionicModelGPUBackendParity/README.md).
It compares every physiological ionic model in scalar, batched-host, and
batched-CUDA execution. The JSON matrix retains all scalar and
`compactBatched` model/tissue rows; CUDA selection is a runtime property of a
batched case with a visible device, not a separate ionic-model name.

The study definition is intentionally limited to case inputs. Run scheduling,
case materialisation, and result storage are supplied externally. Keep all
generated traces and reports in a local ignored results directory.

After two matching host-batched and CUDA-batched cases complete, compare their
exported traces with:

```bash
python3 tutorials/electrophysiologyProtocols/singleCell/setup/studies/ionicModelGPUBackendParity/postprocess_backend_parity.py \
    /path/to/reference-case /path/to/candidate-case
```

The postprocessor checks exported rates plus direct state/current columns at
the saved output precision. Scalar-versus-batched waveform interpretation is
kept separate: distinct ODE integration methods need not be pointwise equal.
