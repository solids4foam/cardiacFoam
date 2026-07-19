# Robustness, Security, and Portability Review

Independent read-only review of commit `f6b798807e8b8c7b685766d6928c16dc5997e150` on branch `no-frontend-minor-errors`. Production code was not modified. The working tree was already dirty, so conclusions below use tracked files and do not treat untracked tutorial material as repository-owned baseline.

## Confirmed findings

### RS-01 — Top-level build masks nonzero `wmake` exits

- **Location:** `Allwmake:4-5`, `Allwmake:37-40`, and the heuristic success check at `Allwmake:42-56`.
- **Evidence:** The script enables `set -e` but not `pipefail`. Both build stages execute `./Allwmake ... | tee log.Allwmake`; Bash therefore uses `tee`'s status for each pipeline. A local semantic reproduction, `bash -c 'set -e; (false 2>&1 | tee /tmp/log); echo masked_rc=$?'`, printed `masked_rc=0`. The later check searches only the literal strings `Error ` and `Stop.`, so a compiler/tool crash, killed process, missing command, or other nonzero exit without either string is reported as success.
- **Concrete failure scenario:** `wmake` is terminated by an out-of-memory condition or a wrapper exits nonzero with `fatal:`/`Killed`; `tee` succeeds, neither grep matches, and `./Allwmake` prints “There were no build errors” and exits zero although a library or solver was not built.
- **Local contract:** `Allwmake:4` explicitly says “Stop at first error”; the root README presents `./Allwmake` as the build entrypoint. Subprocess/build failures must propagate through that entrypoint.
- **Severity / confidence:** **S1 — High / high**.
- **Minimal remediation:** Enable `set -o pipefail` (or use `set -euo pipefail` if unset-variable compatibility is first audited) before the pipelines. Retain log scanning only as supplemental diagnostics, not as the authority for success.
- **Validation:** Run a shell test with a stub `src/Allwmake` that exits a distinctive nonzero status while emitting text containing neither `Error ` nor `Stop.`; assert root `Allwmake` exits nonzero and does not start the applications stage. Then run a normal lightweight and full-mode build.

### RS-02 — Destructive tutorial cleaners act on the caller's directory

- **Location:** The following tracked scripts do not anchor themselves with `cd "${0%/*}"`/`BASH_SOURCE` before invoking `cleanCase` or relative `rm -rf`: `tutorials/NiedererEtAl2011/NiedererEtAl2011verification/Allclean:3-8`; `tutorials/NiedererEtAl2011/electroMechanicalNiedererEtAl2011/Allclean:6-13`; `tutorials/NiedererEtAl2011/purkinjeNiedererEtAl2011/Allclean:3-9`; `tutorials/electrophysiologyProtocols/cableProtocol/eikonal1DCableCV/Allclean:3-8`; `tutorials/electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/Allclean:3-8`; `tutorials/electrophysiologyProtocols/ionicHeterogeneityProbe/Allclean:3-9`; `tutorials/electrophysiologyProtocols/restitutionCurves_s1s2Protocol/Allclean:3-7`; `tutorials/electrophysiologyProtocols/rotorInstability/Allclean:3-9`; `tutorials/electrophysiologyProtocols/singleCell/Allclean:3-10`; `tutorials/manufacturedSolutions/bathBidomain/Allclean:4-11`; `tutorials/manufacturedSolutions/bidomain/Allclean:4-17`; `tutorials/manufacturedSolutions/eikonalECG/Allclean:3-13`; `tutorials/manufacturedSolutions/monodomainNonOrthoMMS/Allclean:4-18`; and `tutorials/manufacturedSolutions/monodomainPseudoECG/Allclean:4-18`.
- **Evidence:** Every removal target in these files is relative to the process CWD. Several use broad patterns such as `processor*`, `[0-9]*`, `0.*`, `constant/polyMesh`, and Bash extglob `!("0"|*[^0-9]*)`. The nearby maintained `tutorials/manufacturedSolutions/monodomainTotalLagrangianEM/Allclean:6` establishes the safe local pattern: `cd "${0%/*}" || exit 1` before cleanup.
- **Concrete failure/abuse scenario:** From an unrelated OpenFOAM case, a user runs `/path/to/cardiacFoam/tutorials/manufacturedSolutions/bidomain/Allclean`. It deletes the unrelated case's time directories, `processor*`, logs, and `constant/polyMesh`. A workflow tool invoking an absolute `Allclean` without setting `cwd` has the same result. This is accidental data loss, not a remote-security threat.
- **Local contract:** Tutorial cleanup must be confined to its own case. The driver normally sets case-local `cwd` (`workflow_runner.py:177-184,223,255-260`), but standalone executable scripts must also be location-independent; the repository already follows this rule in `monodomainTotalLagrangianEM/Allclean` and the root entrypoints.
- **Severity / confidence:** **S1 — High / high**.
- **Minimal remediation:** Add the quoted, fail-closed self-directory `cd` used by `monodomainTotalLagrangianEM/Allclean` to every affected script before sourcing `CleanFunctions` or removing files. Do not broaden or otherwise normalize the deletion patterns in the same patch.
- **Validation:** For each cleaner, copy a minimal case fixture to a path containing spaces; create similarly named sentinel paths in a separate caller CWD; invoke the cleaner by absolute path; assert only the copied case is cleaned and all caller-CWD sentinels remain. Also run each cleaner from its own directory to preserve normal behavior.

### RS-03 — Root build and clean entrypoints break when the checkout path contains whitespace

- **Location:** `Allwmake:2` and `Allwclean:2` use unquoted `cd ${0%/*}`. The correctly quoted local reference is `tutorials/manufacturedSolutions/monodomainTotalLagrangianEM/Allclean:6`; `etc/resolveSolids4Foam.sh:16-19` also consistently quotes derived repository paths.
- **Evidence:** Copying `Allwmake` under a temporary `space dir` and invoking it by absolute path fails at line 2 with `cd: .../space: No such file or directory` and exit status 1. Word splitting occurs before `cd`. The same expression is present in `Allwclean`.
- **Concrete failure scenario:** A user clones the project to `/home/user/Cardiac Research/cardiacFoam`. Both documented root build and clean commands fail before environment checks or any useful work; full and lightweight modes are equally affected.
- **Local contract:** Repository entrypoints deliberately change to their own directory so they are callable from any CWD. That path must be treated as one argument, consistent with the resolver and the newer tutorial cleaner.
- **Severity / confidence:** **S2 — Medium / high**.
- **Minimal remediation:** Change both lines to `cd "${0%/*}" || exit 1` (or the existing `BASH_SOURCE`-based absolute-root idiom if maintainers want symlink behavior standardized).
- **Validation:** Invoke both scripts via absolute and relative paths from outside a checkout whose parent and repository names contain spaces. For the clean test use stubbed OpenFOAM cleanup commands/fixtures so it is non-destructive.

### RS-04 — Top-level clean reports success after subordinate cleanup failures

- **Location:** `Allwclean:20-24` runs library and application cleaners in subshells, but the script has no `set -e` and does not inspect their statuses; `Allwclean:39-41` then ends with `find`, whose success becomes the script's exit status.
- **Evidence:** A failure of `(cd src && ./Allwclean)` or `(cd applications && ./Allwclean)` does not stop execution. Unless the final `find` itself fails, the root command returns zero. This differs from the explicit fail-fast contract used by the build entrypoint at `Allwmake:4-5`.
- **Concrete failure scenario:** `src/Allwclean` fails because `wclean` is missing or a protected build artifact cannot be removed. Root `Allwclean` proceeds, deletes logs, and exits zero, leaving stale objects that can contaminate the next build while automation believes cleanup succeeded.
- **Local contract:** The top-level cleaner orchestrates complete cleanup of libraries and applications (`Allwclean:20-24`); an incomplete clean must be visible to callers. This is especially important when validating full versus lightweight rebuilds.
- **Severity / confidence:** **S2 — Medium / high**.
- **Minimal remediation:** Add fail-fast status propagation (`set -e`, with explicit handling only for intentionally optional cleanup), or collect failures and return nonzero after attempting all independent clean stages. Keep the optional missing `tutorials/Allclean` behavior unchanged.
- **Validation:** In an isolated fixture or with command stubs, make each subordinate `Allwclean` fail in turn and assert root `Allwclean` returns nonzero. Verify a normal full and lightweight clean still succeeds.

## Reviewed areas with no confirmed defect

- Driver workflow subprocesses use argv separation, command allowlisting, case-root-confined resolved CWDs, return-code checks, and timeouts (`core/runtime/workflow.py`, `workflow_runner.py`). The macOS DYLD shell wrapper quotes interpolated values with `shlex.quote`; no command-injection finding was substantiated.
- Run-document execution resolves symlinks and supports an opt-in `DRIVERFOAM_ALLOWED_RUNS_ROOT` boundary (`run_document_exec.py:191-255`), with negative tests for traversal and disallowed commands.
- Sweep expansion caps case counts before materialization and rejects slash/`.`/`..` case identifiers (`sweep_expansion.py:145-156,218-233`). Given the supported Unix/OpenFOAM environments, omission of Windows-specific path rules was not treated as a realistic vulnerability.
- Workflow state and sweep manifests use temporary siblings plus `os.replace` (`workflow_runner.py:66-70`; `sweep_manifest.py:72-76`). Append-only JSONL event/audit files are not process-atomic across multiple concurrent writers, but no repository contract promising concurrent writers was found, so this remains an investigation question rather than a finding.
- Manifests contain absolute host paths by design for local provenance. No secrets or complete environment snapshots are serialized. Their portability tradeoff is visible, but there is no documented relocatable-manifest contract to justify a defect finding.

## Scope note

These are local/operator robustness risks in research tooling. No credible remote attack surface, privilege escalation, secret exposure, or arbitrary command injection was identified in the reviewed project-owned code. The most severe items are data-loss and false-success hazards under plausible developer and automation usage.
