# driverFOAM `--fresh` flag

## Problem

`run --strict` and `sweep-run` silently resume from an existing
`workflow_state.json` / `sweep_manifest.json`. If a case directory from a
previous session is still present and its recorded status is `completed`,
the run reports success and exits 0 without invoking the solver (or, for
`sweep-run`, without even launching the case subprocess) — see
`core/runtime/workflow_orchestrator.py:91` (`while ... status == "pending"`)
and `core/runtime/sweep_runner.py:237` (`elif prior_status == "completed":
outcome = "skipped"`). There is no warning and no flag to force a real rerun.
This was hit in practice on 2026-08-05: a sweep re-run after a solver code
change silently reported yesterday's numbers as fresh results. It was only
caught because the "new" errors matched the old ones to six significant
figures — two different code versions cannot agree that precisely, so
identical numbers meant identical (non-)execution, not agreement.

Manually `rm -rf`-ing the case directory before every re-run works but is
easy to forget and invalidates any before/after comparison silently when
skipped.

## Goal

Add a `--fresh` CLI flag to `action=run`, `action=step`, and
`action=sweep-run` that deletes the resolved output directory before
running, so the workflow executes exactly as it would on a first run. No
change to existing resume behavior when `--fresh` is not passed.

## Semantics

`--fresh` means: delete the resolved output directory tree entirely, then
proceed exactly as the existing "no prior state" code path already does.
This is deliberately whole-tree, not selective — a solver code change
invalidates every artifact under that output directory equally (mesh,
fields, logs, state), not just `workflow_state.json`. Deleting only the
state file was considered and rejected: the orchestrator has no
artifact-existence checks, so every step *would* re-execute, but stale
outputs a differently-scoped rerun doesn't happen to regenerate (leftover
high-timestep directories from a longer previous run, stale postProcessing
files) would silently linger and could be mistaken for current results —
the same class of bug this flag exists to eliminate.

## Hook points

- `run` / `step` (`cli.py`): `_context_from_entry` / `_context_from_run_document`
  resolve `context.output_dir` from the case/tutorial layout (driver-computed,
  not a raw path the user types for deletion). Right after that resolution
  and before `_dispatch_context` is called: if `args.fresh` and the safety
  checks below pass, delete `output_dir` if it exists.
- `sweep-run` (`core/runtime/sweep_runner.py::sweep_run`): right after
  `output_dir = Path(output_dir)` and before `sweep_manifest.json` is read
  (line ~177). Deleting here means the manifest read finds nothing,
  `existing_status_by_case` stays empty, and every case naturally takes the
  existing full-fresh-run branch — no change needed to the per-case
  skip/retry logic itself.

## Safety guards (apply before any deletion, both hook points)

`sweep-run`'s `--output-dir` is a raw user-supplied string, unlike
`run`/`step`'s driver-computed path, so these guards matter most there —
applied uniformly for defense in depth:

1. **Path floor**: refuse if the resolved path is the filesystem root, the
   user's home directory, or has fewer than 3 path segments beneath the
   filesystem root (e.g. `/Users/simaocastro` is refused, `/Users/simaocastro/anything`
   is not).
2. **Allowed-root boundary**: if `DRIVERFOAM_ALLOWED_RUNS_ROOT` is set
   (existing env var, already enforced in
   `core/runtime/run_document_exec.py:236-251`), require the resolved path
   be inside it.
3. **Content-signature check** (new, added in response to review): if the
   target directory exists, refuse deletion unless it contains at least one
   recognizable driverFOAM artifact — `workflow_state.json`,
   `sweep_manifest.json`, or `run_document.json` — within a bounded search
   (top-level plus one level of case subdirectories, i.e. `output_dir/*` and
   `output_dir/*/*`; these markers always live at the case root or the
   sweep's per-case root, so this depth is sufficient). Deliberately not an
   unbounded recursive walk — sweep output dirs can contain hundreds of case
   subdirectories with large mesh trees, and an unbounded scan would be slow
   for no added safety. A directory that doesn't look like driverFOAM's own
   output is never touched, even if it passes guards 1 and 2 — this catches
   a mistyped `--output-dir` pointed at an unrelated directory, which the
   root/depth checks alone would not catch. A directory that doesn't exist
   yet always passes (nothing to lose).
4. **No interactive prompt** (per prior decision — `sweep-run` drives many
   cases unattended), but the exact path about to be deleted is printed
   before deletion, as the only audit trail.

Any guard failure returns the same JSON-error shape the CLI already uses
elsewhere (`status: failed`, `error: <reason>`), exit code 1, no deletion.

## Mutual exclusivity

`--fresh` and `--retry-failed` together is a contradiction (retry-failed
means resume and rerun only failures; fresh means wipe everything and rerun
all). Reject the combination in `_validate_args` (`cli.py:650-681`), same
place `--retry-failed`'s existing `action != "sweep-run"` check lives.

## Usage guidance (for `AGENT_GUIDE.md` and `--help` text)

1. Treat any `--output-dir` passed with `--fresh` as fully disposable —
   never point it at a location that also holds hand-curated files not
   generated through driverFOAM.
2. `--fresh` guarantees deletion, not backup. To keep a before/after
   comparison, copy what you want to keep out of `output_dir` *before*
   rerunning with `--fresh`.
3. For `sweep-run`, the whole `output_dir` tree is wiped, not just cases
   affected by a given change — intentional, since a solver/code change
   invalidates every case's results equally. To preserve some cases
   untouched, use a separate `--output-dir` per logical group rather than
   co-mingling with ones you plan to `--fresh`.
4. Read the printed "about to delete: `<path>`" line the first time you
   point `--fresh` at a new output-dir shape.
5. Any re-run intended as a before/after comparison after a code change
   MUST use `--fresh` (or an equivalent manual `rm -rf`) — otherwise a
   leftover `workflow_state.json`/`sweep_manifest.json` reporting
   `completed` will silently replay the old results.

## Testing

- Reproduce the actual incident at both levels: a case directory with a
  stale `completed` `workflow_state.json` + `--fresh` on `run`/`step` →
  assert the solver command is actually invoked, not skipped.
- Same at the sweep level: a `sweep_manifest.json` claiming all cases
  `completed` + `sweep-run --fresh` → assert every case is re-materialized,
  re-planned, and re-run.
- Guard-refusal cases: home directory, filesystem root, shallow-depth path,
  outside `DRIVERFOAM_ALLOWED_RUNS_ROOT` when set, and a directory that
  exists but contains no driverFOAM artifact markers.
- `--fresh` + `--retry-failed` rejected by argument validation.
- Regression: existing resume behavior is unchanged when `--fresh` is not
  passed.

## Documentation

Update `AGENT_GUIDE.md` with the resume trap and `--fresh` as the fix,
alongside the existing `wclean` build-staleness note.
