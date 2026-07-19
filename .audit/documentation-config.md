# Agent 9 — Documentation and Configuration Audit

## Scope and baseline

- Role: Documentation and Configuration Advocate (independent discovery).
- Repository/branch: `cardiacFoam`, `no-frontend-minor-errors`.
- Commit reviewed: `f6b798807e8b8c7b685766d6928c16dc5997e150`.
- Review basis: committed (`HEAD`) root/component READMEs, tutorial instructions,
  utility help declarations, runtime parsing, build manifests, and build-mode
  resolver logic. Untracked documentation and locally modified production files
  were excluded from conclusions.
- File classification: all files cited below are hand-maintained project code or
  hand-maintained project documentation; none is generated, vendored, or in the
  `solids4foam` submodule.
- Production code was not edited.

## Findings

### DOC-001 — `checkMeshGeometry` documentation contradicts the executable CLI and safety default

- **Severity:** S2 — Medium
- **Confidence:** high
- **Documentation:** `applications/utilities/checkMeshGeometry/README.md:3-15,19-26,31-33`
- **Executable contract:** `applications/utilities/checkMeshGeometry/checkMeshGeometry.C:50-67,95-125,132-150`
- **Evidence:**
  - The README says the utility “rescales it in-place,” says an automatically
    detected mismatch overwrites `constant/polyMesh/points`, and presents plain
    `checkMeshGeometry` as the writing invocation.
  - The executable deliberately makes detect-only behavior the default:
    `doWrite` is true only when `-scale` or `-rescale` is present
    (`checkMeshGeometry.C:123-125`).
  - The README advertises `-noScale`, but the executable registers no such
    option. It registers `-rescale`, `-region`, and `-scale`
    (`checkMeshGeometry.C:50-67`).
  - The README classifies dimensions from `1` through `999` as millimetres. The
    executable's maintained threshold is `20 <= maxDim < 1000`, explicitly to
    avoid mistaking whole-torso SI meshes for millimetres
    (`checkMeshGeometry.C:95-111`).
  - The README says only the default mesh region is read, while the executable
    supports `-region` (`checkMeshGeometry.C:56-61,74-87`).
- **Violated local contract/canonical reference:** OpenFOAM CLI help is defined by
  the executable's `argList::add*Option` calls, and the implementation comments
  explicitly establish opt-in writing as the safety policy. A utility README
  must reproduce that public interface rather than invert it.
- **Impact:** A maintainer following the README can run plain
  `checkMeshGeometry` expecting a mesh correction but receive only a warning,
  then try an invalid `-noScale` option. The incorrect 1 m threshold also tells
  users that valid metre-scale torso meshes are millimetre meshes. This can
  block setup or encourage an erroneous manual scale operation.
- **Minimal remediation:** Rewrite the summary/table/usage/options to state that
  detection is the default; document `-rescale`, `-scale <factor>`, and
  `-region <name>`; remove `-noScale`; and set the documented automatic bands to
  `<20` (metres), `[20,1000)` (millimetres), and `[1000,1e6)` (micrometres),
  with the executable's behavior at or above `1e6` described accurately.
- **Required validation:**
  1. Build the utility and compare `checkMeshGeometry -help` with the README.
  2. On disposable meshes with maximum dimensions 1, 20, 999, 1000, and 1e6,
     verify detection at boundary values.
  3. Hash `constant/polyMesh/points` before/after plain invocation (unchanged),
     `-rescale` (auto-scaled), and `-scale` (explicitly scaled).

### DOC-002 — `runPurkinjeGraph` advertises nonexistent `-nSteps` and `-deltaT` options

- **Severity:** S2 — Medium
- **Confidence:** high
- **Documentation:** `applications/utilities/runPurkinjeGraph/README.md:11-15,23-29`
- **Executable contract:** `applications/utilities/runPurkinjeGraph/runPurkinjeGraph.C:99-118,144-179`
- **Evidence:** The README says execution advances for `nSteps`, lists
  `-nSteps <N>` with a default of 10000, and lists `-deltaT <dt>`. The executable
  registers only the utility-specific `-conductionDomain` option
  (`runPurkinjeGraph.C:107-113`). Its loop uses `runTime.run()`,
  `runTime.deltaTValue()`, and `runTime.endTime()`, so duration and step size are
  controlled by `system/controlDict`, not by those documented CLI switches
  (`runPurkinjeGraph.C:144-179`).
- **Violated local contract/canonical reference:** The executable's `argList`
  registrations and `Time` loop are the canonical public CLI/configuration
  contract. The nearby `-conductionDomain` row correctly mirrors an actually
  registered option and is the canonical documentation pattern.
- **Impact:** A diagnostic run scripted from the README fails during argument
  parsing instead of running. Even if a user omits the invalid flags, the stated
  10000-step default creates a false expectation about runtime and output size.
- **Minimal remediation:** Remove `-nSteps` and `-deltaT` from the options table;
  state that `startTime`, `endTime`, `deltaT`, and write scheduling come from
  `system/controlDict`; revise step 3 accordingly. Do not add CLI aliases merely
  to preserve never-implemented documentation unless maintainers separately
  approve that public feature.
- **Required validation:** Build the utility, compare `runPurkinjeGraph -help`
  with the revised options table, and run a minimal graph case with distinctive
  `controlDict` `deltaT`/`endTime` values to confirm the reported values and step
  count.

### DOC-003 — Committed tutorial indexes and the Purkinje tutorial README name paths that do not exist

- **Severity:** S3 — Low
- **Confidence:** high
- **Documentation:** `README.md:51-60`; `tutorials/README.md:16-22`;
  `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/README.md:1,10-28`
- **Executable/repository reference:** `tutorials/Alltest-regression:23-34,45-61`;
  committed directory `tutorials/electrophysiologyProtocols/`; committed
  directory `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/`
- **Evidence:**
  - The root index names a tutorial group `singleCellprotocols/`, but no such
    committed directory exists. The cases it describes live under
    `electrophysiologyProtocols/`, as also shown by the regression driver
    (`tutorials/Alltest-regression:25,33`).
  - `tutorials/README.md:18` names
    `NiedererEtAl2011/monodomainPurkinjeNiedererEtAl2011`, but the committed and
    regression-tested path is
    `NiedererEtAl2011/purkinjeNiedererEtAl2011`
    (`tutorials/Alltest-regression:28`).
  - The tutorial's own README repeats the obsolete longer name in its title and
    folder tree (`.../purkinjeNiedererEtAl2011/README.md:1,13`).
- **Violated local contract/canonical reference:** Committed directory names and
  the executable regression driver's `REGRESSION_TESTS` paths establish the
  runnable tutorial layout. Documentation intended as a directory index must
  use those exact paths.
- **Impact:** Copy/paste navigation and scripted `cd` commands based on the
  indexes fail, and a new maintainer cannot reliably map the documented
  Purkinje tutorial to the case exercised by regression testing.
- **Minimal remediation:** Rename the root table entry to
  `electrophysiologyProtocols/`; replace both occurrences of
  `monodomainPurkinjeNiedererEtAl2011` with `purkinjeNiedererEtAl2011`. If the
  longer title is retained as a conceptual case name, clearly distinguish it
  from the filesystem path.
- **Required validation:** Check every backticked tutorial path in the two index
  READMEs against `git ls-tree -d`; run a link/path checker in CI; confirm
  `tutorials/Alltest-regression` copies the corrected Purkinje path.

## Reviewed areas with no reportable discrepancy found

- The scalar, batched, and manufactured model lists in
  `src/ionicModels/README.md:114-152` match
  `src/ionicModels/Make/files:6-35` at the reviewed commit.
- The top-level `myocardiumSolver` names in `src/electroModels/README.md:22-30`
  agree with the current layered selection architecture.
- The full/lightweight compilation split described at a high level in the root
  and electromechanical READMEs agrees with `src/Allwmake:22-31` and
  `etc/resolveSolids4Foam.sh:27-95`. The root documentation is terse, but no
  executable contradiction was established, so this audit does not elevate
  terseness alone to a finding.

## Discovery validation performed

- Inspected committed content with `git show HEAD:<path>` and committed layout
  with `git ls-tree`; this avoided treating unrelated working-tree files as
  repository documentation.
- Compared documented options to `argList::addOption`, `addBoolOption`, and
  runtime consumption sites.
- Compared tutorial path claims to `tutorials/Alltest-regression` and committed
  directories.
- Compared the documented ionic-model inventory to `src/ionicModels/Make/files`.
- No builds or numerical tutorial runs were needed to establish these static
  documentation/CLI discrepancies; the runtime checks listed per finding remain
  required after repair.
