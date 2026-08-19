# Publication Readiness Review

Review date: 2026-08-15. This assessment refers to the inspected working tree,
not only to committed code.

## Verdict

driverFOAM is close to supporting a credible research-software paper if it is
presented as a cardiacFoam/OpenFOAM orchestration system. It is not yet ready for
a strong scientific-methods paper or a broad autonomous-agent claim.

The publishable contribution is not a replacement for the C++ solver. The same
C++ executable performs the numerical computation. The prospective contribution
is a typed and inspectable layer for:

- case and configuration resolution;
- encoded pre-flight consistency checks;
- workflow-DAG construction and command authorization;
- parameter sweeps and execution state;
- failure localization and recovery;
- provenance and artifact contracts;
- deterministic tools through which a human, CI system, or language agent can
  operate OpenFOAM.

Consequently, driverFOAM should be expected to reproduce the direct solver's
numerical results. Its benefits must be established using operational metrics:
invalid runs prevented, diagnosis and setup time, recovery success, provenance
completeness, reproducibility, and orchestration overhead.

## Current evidence

### Strengths

- The local suite collected 1,333 tests across 105 test modules.
- In the inspected environment, 1,321 tests passed and 11 skipped.
- There are 16 declared verification experiments, nine canonical regression
  cases, C++/Python catalog drift guards, a RunDocument model, workflow security
  checks, provenance checks, and standalone-package CI.
- The cardiacFoam build workflow exercises OpenFOAM v2312, v2412, and v2512 in
  two build configurations.
- Manufactured-solution and benchmark infrastructure already exists for
  monodomain, bidomain, eikonal/ECG, bath-bidomain, Purkinje/1D-3D coupling, and
  the Niederer slab.
- Security limitations are documented rather than hidden: plugins, case
  scripts, dictionary code directives, binaries, and the local environment are
  trusted rather than sandboxed.

### Blocking evidence gaps

1. **The test tree is not fully green.** One tutorial-characterization drift
   guard failed. The current digests differ from the committed fixture for
   several manufactured tutorials. The changes may be intentional, but they
   require review and an explicit fixture update or code correction.
2. **Solver-dependent driver equivalence is not a routine publication gate.**
   The repository has a phase-2 agent/driver regression harness, but it
   self-skips without a built and sourced solver and is not visibly run by the
   normal driverFOAM CI job.
3. **Experiment declarations exceed archived results.** Seven of the 16
   experiments use `driver_sweep`, while nine still use legacy Bash. Thirteen
   declare reference CSVs, but only one expected result CSV was present in the
   inspected working tree.
4. **One CI path appears stale.** `.github/workflows/buildAndTest.yml` refers to
   `openfoam_driver/tests/test_apply_overrides.py`, while the test is currently
   under `openfoam_driver/tests/core/`.
5. **The full data model is not yet project-neutral.** RunDocument v2 and its
   validation vocabulary remain cardiac-shaped. An external plugin can exercise
   the generic execution mechanisms, but it cannot yet define a fully native
   project configuration schema.
6. **No out-of-tree external project proves the plugin claim.** Structural
   protocols and minimal tests are useful, but they are not empirical evidence
   that another OpenFOAM project integrates without core changes.
7. **No controlled user or agent study establishes improved task performance.**
   Architectural plausibility and unit tests cannot substitute for measured
   success rate, effort, cost, safety, or time.

## Defensible and indefensible claims

| Topic | Defensible now | Not defensible yet |
|---|---|---|
| Solver behavior | driverFOAM orchestrates existing OpenFOAM executables | driverFOAM improves numerical accuracy or solver speed |
| Validation | strict planning detects a declared set of structural and semantic inconsistencies | strict planning guarantees physical correctness |
| Generality | major execution mechanisms are solver-neutral within OpenFOAM | the complete current package is project- or simulator-agnostic |
| Security | RunDocument ingestion and workflow commands have explicit checks | plugins, cases, or dictionary contents are sandboxed |
| Reproducibility | selected inputs, plugin identity, workflow state, and artifacts are recorded | every dependency and environment detail needed for independent reproduction is captured |
| Agents | driverFOAM exposes deterministic, machine-readable tools suitable for agents | an autonomous agent performs reliable scientific discovery |

## Recommended research question

> Can a typed, plugin-mediated planning and execution layer reduce invalid
> OpenFOAM runs and improve reproducibility, recovery, and parameter-study
> automation for cardiac simulations without changing numerical solutions or
> imposing material runtime overhead?

### Hypotheses

- **H1 — Numerical non-interference:** identical inputs executed directly and
  through driverFOAM yield equivalent outputs within predeclared tolerances.
- **H2 — Defect detection:** driverFOAM prevents and correctly localizes more
  invalid launches than direct OpenFOAM, shell, and minimal Python baselines.
- **H3 — Recovery:** driverFOAM resumes interrupted studies with less repeated
  computation and no stale-result acceptance.
- **H4 — Productivity:** users or agents complete setup, sweep, diagnosis, and
  replay tasks with fewer interventions and less time.
- **H5 — Overhead:** orchestration overhead is negligible relative to realistic
  solver workloads and scales predictably with sweep size.
- **H6 — OpenFOAM project generality:** an independently packaged solids4foam
  plugin operates against an unmodified released driverFOAM core.

## Baseline design

“C++ alone” is too ambiguous for the paper. Use four explicit baselines:

1. **Direct OpenFOAM:** hand-edited case plus its canonical `Allrun` or direct
   solver command.
2. **Minimal Python:** a small `subprocess`-based runner and parameter loop with
   no driverFOAM services.
3. **driverFOAM deterministic:** the planner/executor without any LLM.
4. **Agent conditions, only if claimed:** the same language model with raw
   shell/OpenFOAM tools and with driverFOAM tools.

The deterministic driver condition isolates the software contribution. The two
agent conditions isolate whether the interface actually helps an agent.

## Why this evaluation matches publication practice

Workflow papers generally support generality using reproducible execution,
heterogeneous platforms, fault handling, provenance, scaling, and real use
cases—not by requiring every unrelated simulator. Snakemake demonstrated a
domain workflow that moved from workstation to cluster [R1]. Nextflow focused
on portable and reproducible execution [R2]. Pegasus supported broader claims
with multiple scientific domains, distributed platforms, large workflows,
failure handling, and provenance [R3]. A controlled workflow-system comparison
used expressiveness, modularity, scalability, robustness, reproducibility,
interoperability, and ease of development as evaluation dimensions [R4]. signac
demonstrated lightweight workflow/data management and research use cases [R5,
R6].

For an agent claim, task-level executable evaluation is necessary.
ScienceAgentBench uses expert-validated tasks, multiple models and scaffolds,
repeated attempts, execution/results metrics, and cost [R8]. MetaOpenFOAM and
related OpenFOAM agent work provide closer CFD precedents, but some remain
preprints and should be presented as such [R9–R11].

## Recommended paper type and scope

### Near-term: research-software paper

Recommended title:

> **driverFOAM: Typed Planning, Validation, and Reproducible Workflow
> Orchestration for Cardiac OpenFOAM Simulations**

This paper should emphasize software architecture, numerical non-interference,
defect detection, recovery, provenance, and one independently reproducible case.

### Stronger methods/workflow paper

This requires the complete controlled benchmark: frozen fault corpus, direct
and Python baselines, quantitative recovery/overhead results, an independent
reproduction, and external-plugin evidence.

### Agent paper

Make this a separate paper or a clearly isolated experimental section. It needs
a hidden task set, multiple model families and seeds, raw-shell and deterministic
baselines, component ablations, unsafe-action reporting, expert scoring, and
cost/time results. Do not let the stochastic agent obscure the deterministic
driverFOAM contribution.

## Suggested manuscript structure

1. Problem statement and claim boundaries.
2. Related work: OpenFOAM automation, scientific workflow systems, research
   software reproducibility, and agents for science.
3. Architecture and trust model.
4. Research questions, hypotheses, cases, baselines, metrics, and statistical
   analysis plan.
5. Numerical non-interference and manufactured-solution results.
6. Fault-injection detection and localization results.
7. Recovery, provenance, and replay results.
8. Sweep overhead and scaling results.
9. solids4foam external-plugin case study.
10. Optional human/agent task study.
11. Limitations and threats to validity.
12. Artifact and data availability.

## Venue fit

- **OpenFOAM Journal:** strongest domain fit when the paper emphasizes concrete
  OpenFOAM benefit, executable cases, and numerical evidence [R17, R18].
- **Journal of Open Source Software:** suitable for a compact mature-software
  paper once public development history, research use, release, archival DOI,
  tests, CI, documentation, and contribution guidance meet its requirements
  [R19]. The 2025 solids4foam JOSS paper is a directly relevant precedent [R12].
- **SoftwareX or Journal of Open Research Software:** appropriate for a
  reuse-focused software description with public artifacts [R20, R21].
- **Workflow/computational-science journal:** appropriate only after controlled
  scalability, fault, multi-platform, and reproducibility results.

## Decision rule

The project is ready for a software-paper submission when all release and CI
gates are green, at least one complete driver-driven numerical study is archived
and independently reproduced, and the manuscript's claims match the current
cardiac/OpenFOAM scope.

It is ready for a scientific workflow/methods submission only when the
predeclared experiments in `EXPERIMENTAL_EVIDENCE_TABLE.md` produce quantitative
evidence against explicit baselines.

References [R1–R22] are listed in [`REFERENCES.md`](REFERENCES.md).
