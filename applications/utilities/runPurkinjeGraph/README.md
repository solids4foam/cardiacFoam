# runPurkinjeGraph

Diagnostic utility for running a `conductionSystemDomain` (Purkinje graph)
in isolation, without constructing the myocardium mesh or any coupling.

Use this to check that a Purkinje network evolves correctly before integrating
it into a full coupled run.

## What it does

1. Reads `constant/electroProperties` and selects the `conductionNetworkDomains`
   block from the active solver coefficients.
2. Constructs a `conductionSystemDomain` from the chosen graph dictionary.
3. Advances the domain for `nSteps` time steps at step size `deltaT`.
4. Writes domain output at each `writeControl` interval.

## Usage

```bash
runPurkinjeGraph [options]
```

## Options

| Option | Default | Description |
|---|---|---|
| `-conductionDomain <name>` | first entry | Name of the `conductionNetworkDomains` sub-dictionary to run |
| `-nSteps <N>` | 10000 | Number of time steps to advance |
| `-deltaT <dt>` | from `controlDict` | Override the time-step size (s) |

## Notes

- Requires a valid mesh in `constant/polyMesh` (used by the domain constructor)
  even though no spatial PDE is solved.
- Runs in serial only (`noParallel`).
- No myocardium coupling currents are injected; terminal PVJ buffers are zero
  throughout.
- Output fields follow the domain's normal write rules (`conductionSystemDomain`
  writes nodal `Vm1D`, `activationTime`, and Purkinje edge data).
