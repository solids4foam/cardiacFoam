# CLAUDE.md

Project overview and architecture: [README.md](README.md).

## How cases are built, run, and tested here

Each tutorial under `tutorials/` is self-contained: its own `Allrun` builds
and launches the case, its own `Allclean` resets it, and its own
`regression/regressionTest.sh` checks the result against the tutorial's
reference data. The suite-wide driver is `tutorials/Alltest-regression`,
which sweeps every tutorial's regression check. It requires the
`CARDIAC_REGRESSION_BUILD_MODE` environment variable to be set to exactly
`with-solids4foam` or `lightweight` (it exits non-zero otherwise), matching
how the case was built.

**Do not write a bespoke shell script that sources OpenFOAM's `bashrc`/
`RunFunctions` and mutates tracked case dictionaries directly.** This has
caused real, silent damage before: an ad-hoc script once flipped a tracked
`fvSchemes` default and overwrote `box.geo.template` with the wrong
mesh-generation variant, undetected until a manual audit. If an existing
tutorial script already does what you need, reuse it rather than adding a
new one.

A tutorial's own committed `Allrun` and `Allclean` are shell, and they are
exactly how a case is meant to be built and run. Raw shell is likewise fine
for genuinely one-off, non-simulation glue — the kind already established at
the repo root (`reproduce_verification.sh`, `Allwmake`). What to avoid is a
*new* ad-hoc script that reaches into a case and rewrites its dictionaries.

An optional external orchestration add-on for building, sweeping, and
validating cases exists outside this repository and ships its own
documentation; it is not part of this codebase.
