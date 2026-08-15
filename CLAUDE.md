# CLAUDE.md

Project overview and architecture: [README.md](README.md).

## Always use driverFOAM for cases, sweeps, and runs

`applications/scripts/driverFoam/` (the `foamctl` CLI / `openfoam_driver`
package) is the only supported way to build, validate, launch, and sweep
OpenFOAM cases in this repo. It owns dictionary mutation, the strict
pre-flight planner, artifact manifests, and sweep orchestration
(`foamctl plan --strict`, `foamctl run --strict`,
`foamctl sweep-run --spec sweep.json`). Read
[`applications/scripts/driverFoam/AGENT_GUIDE.md`](applications/scripts/driverFoam/AGENT_GUIDE.md)
before driving any case. A packaged walkthrough of the same workflow is also
available as the `driverfoam-assistant` skill.

**Do not write a bespoke shell script that sources OpenFOAM's `bashrc`/
`RunFunctions` and calls the solver binary directly.** That bypasses
driverFOAM's validated dict-mutation path and re-implements logic it already
provides (mesh setup, parameter sweeps, convergence ladders, artifact
tracking) in untracked bash. This has caused real, silent damage before: an
ad-hoc script once flipped a tracked `fvSchemes` default and overwrote
`box.geo.template` with the wrong mesh-generation variant, undetected until a
manual audit. If an existing tutorial script already does
what you need, reuse it. Otherwise, write a `sweep.json` (or `RunDocument`)
and drive it through `foamctl` — don't add a new `run_*.sh`.

Raw shell is still fine for genuinely one-off, non-simulation glue — the kind
already established at the repo root (`reproduce_verification.sh`,
`Allwmake`) — but not for anything that builds, mutates, or runs a case.
