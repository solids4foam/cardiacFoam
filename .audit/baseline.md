# Audit baseline

- Repository URL: `git@github.com:solids4foam/cardiacFoam.git` (`github` remote)
- Branch: `no-frontend-minor-errors`
- Commit: `f6b798807e8b8c7b685766d6928c16dc5997e150`
- Submodule: `modules/solids4foam` at `0bd882172db292c29bf41c4233d61cfa5f116168`; the submodule working tree is modified and is treated as external, out-of-scope content.
- Host: macOS; no OpenFOAM environment variables are set.
- Compiler commands available: Apple/system `gcc`, `g++`, and `clang++`.
- OpenFOAM commands unavailable on `PATH`: `wmake`, `foamVersion`.
- Python: `3.14.3` at `/opt/homebrew/bin/python3`; `pytest` is available.
- Other available build tools: `cmake`, `ninja`.
- Working tree: substantially dirty before the audit, including modified tracked files and many untracked research/tutorial files. All pre-existing changes belong to the user and must be preserved. Production findings must distinguish committed baseline behavior from pre-existing working-tree changes where relevant.
- Commands expected to be blocked by environment: `./Allwmake` and OpenFOAM tutorial regressions, because no OpenFOAM environment is loaded and `wmake` is unavailable.
- Safe commands available: repository inspection, Git read-only commands, shell/static checks, and Python tests that do not require OpenFOAM.

## Compatibility-sensitive surfaces

- OpenFOAM runtime-selection strings and dictionary keys, especially `myocardiumSolver`, ionic-model names, and their `*Coeffs` dictionaries.
- `Make/files`, `Make/options`, `Allwmake`, and the full/lightweight resolver contract.
- Python CLI arguments, JSON schemas/manifests, tutorial specifications, and artifact contracts.
- Generated scalar/batched ionic model metadata and CPU/GPU parity.
- Tutorial directory names, dictionaries, regression entrypoints, and documented commands.

## Authority and boundaries

- Authority order: maintained code, build files, component documentation, project memory, audit guidance.
- `modules/solids4foam` is a submodule/external boundary and must not be modified.
- Large ionic equation headers and generated CellML outputs must be classified before proposing edits; generator/template fixes are preferred.
- `.audit/` is the only writable area during discovery.
