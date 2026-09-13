# genericWriter

All the output logic, shared by the ionic, active-tension and electro libraries: writing traces and time series, and reading stimuli.

## What's available

| File | What it does |
|---|---|
| `ionicModelIO` | ionic-model output: traces, field export, selected variables |
| `ionicVariableCompatibility` | matches variable names across models, for export and signal lookup |
| `stimulusIO` | reads and evaluates stimuli |
| `StimulusProtocolPOD.H` | a plain-data copy of a stimulus protocol, for batched and GPU ionic models |
| `activeTensionIO` | active-tension output |
| `ecgModelIO` | ECG output |
| `purkinjeModelIO` | Purkinje-network time series |
| `conductivityFieldIO` | resolves the conductivity used by the tissue solvers |

## Folders

All files sit directly in `src/genericWriter/`.

## What this does not own

- The physics models themselves. This library only reads inputs and writes outputs for them.
