# HeartSimTemplate Status

## Done

- `constant/electroProperties` has been aligned to the current template-truth dictionary family:
  - `myocardiumSolver`
  - `monodomainSolverCoeffs`
  - `conductionNetworkDomains`
  - `domainCouplings`
  - `ecgDomains` now placed inside `monodomainSolverCoeffs`
- The nested `driverFOAM` key schema is aligned with the same contract.
- Tutorial contracts are exposed in machine-readable form through driver introspection.
- `HeartSimTemplate` now provides a machine-readable workflow contract:
  - [workflow_contract.json](./workflow_contract.json)
- The current major heart workflow is represented as:
  - `monodomain_purkinje_pseudo_ecg`
- The contract now distinguishes:
  - core required files
  - solver required files
  - conditional files
- The current template is explicitly marked as symbolic authoring input, not a directly runnable case.

## Missing For Full Machine-Ready Use

- Replace symbolic ionic-model placeholders with concrete model names:
  - `<myocardiumIonicModel>`
  - `<purkinjeIonicModel>`
- Replace symbolic ionic export placeholders with real export-variable lists.
- Define a code-backed protocol for Purkinje ionic export discovery, analogous to the myocardium export-variable flow.
- Decide whether `HeartSimTemplate` should remain a generic template folder or become a curated registered tutorial in `driverFOAM`.
- Ground the template against a real full-heart production case/dictionary so the workflow contract is based on a final reference case, not only symbolic templates.

## Current Interpretation For Agents

- Use `dict_entries.py` for key-level schema.
- Use driver `tutorial_contract` for case file inventory and case-parameter metadata.
- Use `workflow_contract.json` for workflow assembly rules and compatibility guidance.
- Treat `HeartSimTemplate` as an authoring template that still requires substitution before execution.
