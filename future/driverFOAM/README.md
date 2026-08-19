# driverFOAM Publication Planning

This folder contains a skeptical, evidence-driven plan for publishing
driverFOAM. It is planning material, not evidence that the proposed experiments
have already passed.

## Documents

- [`PUBLICATION_READINESS.md`](PUBLICATION_READINESS.md): current scientific
  and software-publication verdict, defensible claims, paper questions, and
  recommended paper structure.
- [`EXPERIMENTAL_EVIDENCE_TABLE.md`](EXPERIMENTAL_EVIDENCE_TABLE.md): proposed
  experiments, baselines, metrics, acceptance criteria, and artifacts.
- [`SOLIDS4FOAM_CASE_STUDY.md`](SOLIDS4FOAM_CASE_STUDY.md): design for an
  independently packaged solids4foam plugin and external-validity study.
- [`ROADMAP.md`](ROADMAP.md): ordered implementation and experimental roadmap
  with dependencies and completion gates.
- [`REFERENCES.md`](REFERENCES.md): annotated primary literature and official
  publication/venue guidance.
- [`references.bib`](references.bib): initial BibTeX database for the paper.

## Recommended paper scope

The most defensible near-term scope is:

> driverFOAM is a typed planning, validation, provenance, and workflow
> orchestration layer for reproducible cardiacFoam/OpenFOAM simulation studies.

An OpenFOAM-wide claim should be made only after an out-of-tree solids4foam
plugin passes planning, execution, sweep, resume, artifact, and provenance tests
without changes to driverFOAM core.

## Evidence language

Use these distinctions consistently:

- **Plan-valid:** the case satisfies the structural and semantic checks encoded
  by driverFOAM and its selected plugin.
- **Execution-ready:** required runtime environment and executables are
  available in addition to plan validity.
- **Numerically verified:** results satisfy a separately defined analytical,
  manufactured-solution, benchmark, or regression acceptance criterion.
- **Scientifically validated:** results agree with an appropriate external
  physical or experimental reference for the intended use.

Neither plan validity nor execution readiness proves physical correctness.

