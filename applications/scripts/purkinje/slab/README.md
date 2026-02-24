# Purkinje Slab

This module is the Purkinje slab product-facing entrypoint.

Engine implementation lives in:

- `cardiac_preproc/src/cardiac_preproc/steps/purkinje_slab.py`

Product CLI:

- `slab.py`

Config (colocated with CLI):

- `config_slab.py`

Compatibility wrapper:

- `run_slab.py`

Usage:

```bash
python3 purkinje/slab/slab.py --help
```

Core implementation remains in `cardiac_preproc/src/cardiac_preproc/steps/purkinje_slab.py`.
