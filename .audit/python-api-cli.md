# Python API and CLI discovery audit

## Scope and ownership

Read-only review of project-owned Python under `applications/scripts/driverFoam/`,
`applications/scripts/cellML2foam/`, and representative standalone project scripts.
The `presentation/vendor/` tree was classified as vendored and excluded. Generated
`applications/scripts/driverFoam/cardiacfoam_tutorials_driver.egg-info/` metadata was
used only as corroborating evidence, not treated as a repair target. Tutorial setup
scripts are reference/workflow inputs and were sampled rather than normalized.

The maintained driver package establishes the local conventions: `pathlib.Path` at
filesystem boundaries, dataclasses for stable runtime records, JSON-returning
`to_json()` methods, argument-vector subprocess calls, explicit non-zero CLI returns,
and focused contract tests. The older `cellML2foam` utility is an un-packaged script
suite with implicit current-working-directory contracts and currently has no tests.

## Correctness defects

### PY-01 — `cellML2foam --outdir` cannot receive generated headers

- **Evidence:** `applications/scripts/cellML2foam/cellML2Foam.py:125-129` exposes
  `--outdir`; `cellML2Foam.py:152-158` calls `run_pipeline()` without passing it.
  `applications/scripts/cellML2foam/src/pipeline.py:214-220` consequently passes the
  relative filename `<Model>_<Year>.H` to `run_mapping()`. Finally,
  `applications/scripts/cellML2foam/src/sort_folder.py:85-98` requires both generated
  headers under `outdir` and raises if they are absent. A focused mock invocation of
  `run_pipeline(..., start="ansic", end="openfoam", model="Foo_2020")` observed the
  mapping call `(.../sim.c, "Foo_2020.H")`, confirming that no destination is carried.
- **Local contract/reference:** argparse help at `cellML2Foam.py:125-129` calls this the
  output directory; `sort_folder.py:85-98` consistently interprets it as the generated
  files' source directory as well as the model-folder destination.
- **Impact:** every successful OpenFOAM conversion using a non-current `--outdir`
  writes headers in the process CWD, then terminates with `Expected file not found`,
  leaving partial output in two locations. This makes a public CLI option unusable.
- **Severity/confidence:** **S1 (High), high confidence**.
- **Minimal repair:** add an output-directory/path parameter to `run_pipeline`, create
  it as needed, and pass `outdir / f"{name}_{year}.H"` to `run_mapping`; keep the two
  generated headers and subsequent `sort_folder` operation under that directory.
- **Validation:** add a temporary-directory CLI/unit test with CWD distinct from
  `--outdir`; stub transformation/Myokit input as needed; assert exit 0, both generated
  headers and wrappers under `outdir/<Model>/`, and no generated headers in CWD.

### PY-02 — `run_mapping()` splits paired outputs across directories

- **Evidence:** `applications/scripts/cellML2foam/src/mapping_engine.py:376-381`
  accepts an arbitrary `output_c` path but reduces the names-header output to the bare
  basename. It writes `output_c` at lines 416-424, while line 426 calls
  `generate_header(output_h, ...)` in the current directory. The generated source
  includes that bare sibling name at line 418.
- **Local contract/reference:** `sort_folder.py:85-98` treats `<Model>_<Year>.H` and
  `<Model>_<Year>Names.H` as a co-located pair. The function's return at lines 431-435
  also presents both as outputs from one operation.
- **Impact:** programmatic callers that provide a directory in `output_c` receive one
  file there and the required included header in an unrelated CWD. Compilation and
  cleanup become location-dependent. This is also the lower-level cause amplifying
  PY-01.
- **Severity/confidence:** **S2 (Medium), high confidence**.
- **Minimal repair:** derive `output_path = Path(output_c)` and write the names header
  to `output_path.with_name(f"{model_id}Names.H")`; keep only the basename in the C++
  include and return concrete paths consistently.
- **Validation:** unit-test `run_mapping` with an output path in a temporary nested
  directory and assert both files are siblings and the include resolves locally.

## API-contract inconsistencies

### PY-03 — the installed package advertises a dashboard API that is absent

- **Evidence:** `applications/scripts/driverFoam/openfoam_driver/README.md:438-458`
  documents `foamctl dashboard` and
  `python -m openfoam_driver.dashboard.render_case`. The repository contains no
  `openfoam_driver/dashboard/` package, and `openfoam_driver/cli.py` has no dashboard
  action or dispatch. Executing
  `PYTHONPATH=applications/scripts/driverFoam python3 -m openfoam_driver dashboard --help`
  returned argparse exit code 2. `applications/scripts/driverFoam/pyproject.toml:27-37`
  still publishes a `dashboard` dependency extra. Stale generated egg metadata lists
  deleted dashboard modules, corroborating incomplete removal.
- **Local contract/reference:** current CLI actions are documented at
  `openfoam_driver/README.md:66-75` and implemented by the parser/dispatch in
  `openfoam_driver/cli.py`; the package manifest includes only existing package data.
- **Impact:** users can install a substantial optional dependency set and follow the
  committed instructions, but neither advertised entry point can run. Downstream code
  cannot rely on the documented import path.
- **Severity/confidence:** **S2 (Medium), high confidence**.
- **Minimal repair:** maintainer choice is required: restore the dashboard package and
  CLI wiring if it remains supported, or remove the dashboard extra and the obsolete
  README section. Do not retain a dependency-only compatibility extra without a stated
  migration purpose.
- **Validation:** add a packaging/CLI smoke test that installs each advertised extra
  and invokes every README-listed action/import; if removal is chosen, build wheel and
  inspect metadata to ensure the extra and stale modules are absent.

## Maintainability findings

### PY-04 — `cellML2foam` has no regression boundary for its filesystem-heavy CLI

- **Evidence:** all seven maintained Python modules under
  `applications/scripts/cellML2foam/` have no adjacent test directory or test files.
  They coordinate subprocesses and destructive moves (`pipeline.py:75-126`,
  `sort_folder.py:85-98`) through implicit CWD state. PY-01 and PY-02 are straightforward
  path-contract regressions that the current driverFoam suite cannot exercise.
- **Local contract/reference:** `applications/scripts/driverFoam/openfoam_driver/tests/`
  extensively uses temporary directories and mocked subprocesses for the same kinds of
  filesystem and CLI contracts; the audit run completed **712 passed, 3 skipped, 126
  subtests passed**.
- **Impact:** output-location, partial-write, stage-selection, and subprocess-error
  regressions can ship undetected in a code generator whose outputs become compiled
  scientific code.
- **Severity/confidence:** **S2 (Medium), high confidence**.
- **Minimal repair:** introduce a small pytest suite around pipeline validation,
  paired-output placement, non-default `--outdir`, missing tools, and subprocess
  failure propagation. No broad framework or formatting conversion is needed.
- **Validation:** run the new tests from a neutral CWD and include at least one CLI
  subprocess smoke test in addition to unit mocks.

## Style-only observations (not proposed for default repair)

- `applications/scripts/cellML2foam/src/discovery.py:100-101` uses a bare `except` and
  silently discards discovery failures. This differs from the driver's narrower,
  diagnostic-rich error handling. It should be narrowed if this file is touched, but
  the heuristic deliberately returns `None`, so there is insufficient evidence here
  to classify it independently as a behavior defect.
- `applications/scripts/cellML2foam/src/mapping_engine.py:23-40` emits a Python
  `SyntaxWarning` for the invalid `\*` escape in `LICENSE_HEADER`. Use a raw string or
  escape the backslash when nearby functional work touches the file.
- Import ordering, annotations, and whitespace vary within `cellML2foam`; no wholesale
  Black/Ruff conversion is justified by demonstrated behavior.

## Optional tooling recommendations

1. Add a focused `cellML2foam` pytest job using only temporary directories and mocked
   Myokit subprocesses. This directly prevents recurrence of PY-01/PY-02.
2. Add a lightweight distribution smoke test that builds/installs driverFoam and checks
   documented console actions and optional import paths. This directly detects PY-03
   and stale generated package metadata.
3. If a linter is adopted, initially restrict it to definite runtime hazards such as
   bare exceptions and invalid escape sequences in hand-maintained Python; do not make
   repository-wide formatting a prerequisite.

## Validation performed

- `SKIP_ENV_DIAGNOSTICS=1 python3 -m pytest openfoam_driver/tests/ -q` from
  `applications/scripts/driverFoam`: **passed** (`712 passed, 3 skipped, 126 subtests
  passed in 73.42s`).
- Focused mocked `run_pipeline` path probe: confirmed relative generated output.
- `python3 -m openfoam_driver dashboard --help` with repo package on `PYTHONPATH`:
  **failed as expected with exit 2**, confirming the stale documented action.

No production code was modified during discovery.
