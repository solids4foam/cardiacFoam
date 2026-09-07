# Experimental Evidence Table

This is the proposed evidence matrix. Acceptance criteria should be frozen
before generating the final results. Threshold values marked **TBD** require a
short pilot study followed by preregistration; they must not be selected after
seeing the final comparison.

## Core experiments

| ID | Research question | Cases/data | Conditions and baselines | Primary metrics | Proposed acceptance criterion | Required artifacts |
|---|---|---|---|---|---|---|
| E1 | Does driverFOAM alter numerical results? | Single cell, monodomain, bidomain, eikonal/ECG, bath-bidomain, electromechanics or Purkinje | Canonical direct `Allrun`; minimal Python runner; deterministic driverFOAM | Normalized dictionary diff; input hashes; L1/L2/L∞ field errors; activation-time/ECG/displacement differences; convergence order | Direct and driver results agree within case-specific, predeclared numerical tolerances; observed convergence order agrees within **TBD** | Immutable input bundle, environment manifest, raw fields, summary CSV, comparison script |
| E2 | Does strict planning prevent invalid launches? | Frozen corpus of valid cases plus seeded mutations | Direct OpenFOAM/checkMesh; shell/Python checks; driverFOAM | Precision, recall, specificity, F1, false-positive rate, stage localization, invalid launches prevented | High recall on safety/invalid-launch defects without unacceptable false positives; exact target set after pilot | Mutation manifest, ground-truth labels, diagnostic JSON, confusion matrices |
| E3 | Does driverFOAM recover interrupted work correctly? | Multi-step cardiac workflows and sweeps | Clean rerun; minimal checkpoint script; driver resume | Recovery success, repeated solver work, wall time, stale result acceptance, manifest consistency | 100% stale-result rejection for declared perturbations; less repeated work than clean rerun | Kill/fault schedule, state files, provenance snapshots, execution logs |
| E4 | What overhead does orchestration add? | No-op/mock steps; cheap solver case; realistic solver case; sweeps of 1, 10, 100, 1000 where feasible | Shell loop; minimal Python loop; driverFOAM | Planning time, launch overhead, makespan, throughput, peak RSS, state/provenance bytes, solver-time fraction | Median overhead below **TBD** for realistic cases; complexity trend reported rather than hidden | Raw timing JSON, machine details, repetitions, confidence intervals |
| E5 | Can another researcher reproduce a study? | One complete cardiac study, preferably a manufactured/benchmark case and sweep | Documentation-only reconstruction; archived driverFOAM artifact | Setup success, deviations, time to first result, hashes/tolerance agreement, missing provenance | Independent researcher succeeds from a tagged archive without author intervention beyond documented support | Zenodo archive, container/environment recipe, blinded report, checksums |
| E6 | Is the plugin boundary genuinely OpenFOAM-project general? | solids4foam patch test, one nonlinear solid case, one FSI case | Native solids4foam scripts; external driverFOAM plugin | Core files changed, plugin LOC, integration effort, plan/run/sweep/resume pass rate, output equivalence | Zero driverFOAM core edits; all selected workflows pass; outputs match native regression criteria | Separate plugin repository, release, CI logs, conformance report |
| E7 | Does provenance detect relevant changes? | One stable case replayed after controlled changes | Baseline archive; changed input, plugin, workflow, solver binary/library, OpenFOAM version | Change detection rate, false alarms, resume refusal/acceptance, captured dependency coverage | Every declared required change is recorded and invalidates unsafe reuse; unchanged replay remains accepted | Perturbation table, before/after manifests, expected/actual decisions |
| E8 | Does driverFOAM help users? | Setup, sweep, diagnosis, and resume tasks | Direct workflow; driverFOAM; randomized crossover | Completion, human time, interventions, invalid launches, time to first valid result, SUS/NASA-TLX | Directional and statistically supported improvement on preregistered primary endpoint | Protocol, anonymized event logs, analysis notebook, consent/ethics record if required |
| E9 | Does driverFOAM help language agents? | Hidden valid/invalid tasks stratified by physics and operation | Same LLM raw shell; LLM + driverFOAM; driverFOAM deterministic; ablations | Pass@1/pass@k, numerical-validity score, unsafe actions, interventions, tokens, cost, time | Improvement over same-model raw-shell baseline; deterministic condition reported separately | Pinned prompts/models/dates, task corpus, transcripts, expert-blind scores, cost logs |

## Numerical non-interference case matrix

| Case family | Why it is needed | Suggested observables | Direct reference path |
|---|---|---|---|
| Single cell | ODE/configuration and ionic-model selection | Vm and calcium biomarkers, activation/repolarization timing | Existing single-cell regression reference |
| Monodomain manufactured solution | PDE field accuracy and spatial/temporal convergence | Vm, extracellular proxy/ECG, L1/L2/L∞, observed order | Existing monodomain convergence references |
| Bidomain manufactured solution | Coupled elliptic/parabolic configuration | Vm, extracellular potential, error norms, iterations | Existing bidomain references |
| Eikonal/ECG manufactured solution | Different solver family and derived integral artifact | Activation time, pseudo-ECG, bulk/boundary error | Existing eikonal references |
| Bath-bidomain | Multi-region/interface behavior | Heart/bath fields and interface metrics | Existing bath-bidomain references |
| Purkinje or electromechanics | Graph/1D-3D or multiphysics coupling | Coupling error, activation, displacement/active tension | Existing coupling or Niederer-style reference |

Each case must use the same staged input tree for all conditions. If driverFOAM
materializes dictionaries, compare both semantic values and normalized text. Do
not require byte equality when OpenFOAM formatting changes without semantic
change.

## Fault corpus

| Class | Example mutations | Ground-truth outcome |
|---|---|---|
| Required files | Remove `controlDict`, `fvSolution`, plugin configuration, initial field | Block before execution |
| Dictionary syntax/shape | Unbalanced braces, wrong type, missing required block, invalid enum | Block or clearly classify parser limitation |
| Cross-field semantics | Incompatible solver/ionic model/tissue, invalid coupling selection | Block with correct field localization |
| Workflow structure | Missing dependency, cycle, duplicate ID, invalid `cwd` | Block before execution |
| Command boundary | Unknown command, absolute path, case-local binary shadow, unauthorized script | Block according to documented trust policy |
| Environment | Unsourced OpenFOAM, missing executable/library, wrong app path | Mark execution not ready and refuse run |
| Artifacts | Required output absent, wrong producer, stale output, optional artifact absent | Fail only required-artifact cases; avoid false positives |
| Resume/provenance | Changed input, plugin digest, solver binary, workflow, interrupted state | Refuse unsafe reuse or resume from correct step |
| Sweep | Invalid axis, zip length mismatch, path-unsafe case ID, partial completed sweep | Reject invalid spec; recover valid partial sweep |
| Executable dictionary content | `#codeStream`, `#calc`, coded function object | Report trust limitation; do not claim sandboxing |

At least 20–30% of the final mutations should be independently authored and
hidden from the developer implementing the validator. Otherwise the benchmark
measures rule transcription rather than general fault detection.

## Agent benchmark design

Only run E9 if the paper makes a central agent claim.

- Use at least two model families and at least three independent trials/seeds.
- Freeze model identifier, provider version/date, prompts, tools, temperature,
  maximum turns, token budget, and stopping policy.
- Separate tutorial-derived tasks from held-out tasks to reduce contamination.
- Score structural success, executable success, numerical correctness, and
  safety separately.
- Use blind expert review for ambiguous scientific outputs.
- Include component ablations: no strict planner, no provenance/resume, no
  plugin introspection, and no repair loop.
- Report failures and intervention counts, not only the best attempt.

## Statistical reporting

- Repeat timing experiments sufficiently to report medians and bootstrap 95%
  confidence intervals; retain full distributions.
- Use paired comparisons wherever the same case/task is run in each condition.
- Predeclare one primary endpoint per hypothesis and treat the rest as
  secondary/exploratory.
- Report absolute effects as well as percentages.
- Do not pool novice and expert user results without testing the expertise
  interaction.
- Publish raw machine-readable results and the analysis code used for tables and
  figures.

## Minimum evidence packages

| Submission target | Minimum completed experiments |
|---|---|
| JOSS/JORS software paper | E1, E5, green CI/release evidence; E2 or E6 strongly recommended |
| OpenFOAM Journal full paper | E1, E2, E3, E4, plus domain verification and reproducible artifacts |
| SoftwareX | E1, E2, E4, E5, E6, with reuse documentation |
| Agent/AI-for-science paper | E1–E4 and E9; E9 must include baselines and ablations |

