# driverFOAM workflow definitions: complete characterization

**Date:** 2026-08-20
**Status:** Characterization only. No code changes proposed here, by request.
**Scope:** how a case's execution steps are defined *within driverFOAM*, which
definition runs, and where they disagree.

Every claim was produced by resolving real entries and parsing real files, not
by reading code and inferring. Where a measurement was imprecise, that is said.

---

## 0. The layering, stated first because it bounds everything else

`Allrun` is **OpenFOAM's** convention: the basic run script for a case.
driverFOAM is built on top of it and does not own it, author it, or reconcile
against it.

That means `Allrun` is **not** a competing source of truth. driverFOAM's own
definitions are two:

| definition | form | where | how many cases |
|---|---|---|---|
| factory `workflow_dag` | Python, built by `make_spec` | `plugins/cardiacfoam/tutorials/*.py` | 13 registered |
| `workflow_contract.json` | declarative JSON `steps` | case root | 6 |

An earlier draft of this document treated `Allrun` as a third definition and
counted "7 divergent tutorials" by comparing DAGs against it. That count is
withdrawn: a DAG differing from a shell script is expected when the DAG is a
layer above it, not a defect.

`Allrun` still appears below, in one legitimate role: as **evidence of the
case author's intent**. When a DAG and a contract disagree about whether a
case needs `blockMesh`, what the case's own run script does is admissible
corroboration. It is evidence, not a specification.

### 0.1 driverFOAM's three legitimate contacts with `Allrun`

1. **Detects** it. `registry._is_case_directory` and
   `is_runnable_without_workflow` treat its presence as evidence that a folder
   is a runnable case.
2. **Invokes** it as a DAG step. `cable1DCVConvergence` has a DAG of exactly
   `[Allrun]`; `cable1DRestitution` has `[Allclean, Allrun, Allrun.post]`.
   The DAG delegates wholesale to the OpenFOAM entry point.
3. **Never modifies** it — with one exception, below.

### 0.2 The one place that violates the layering

`plugins/cardiacfoam/sweep.py` **authors** an `Allrun` when materializing a
sweep case:

```python
allrun_body = "blockMesh\ncardiacFoam\n" if result.get("needs_block_mesh") else "cardiacFoam\n"
allrun_path.write_text("#!/bin/sh\n" + allrun_body)
```

Two things about it are worth noting. It is driverFOAM writing a file it does
not own, contradicting the stated layering. And it is a *raw* script — bare
commands under `#!/bin/sh`, not `RunFunctions` / `runApplication` — so a
materialized sweep case carries an `Allrun` that does not follow the OpenFOAM
convention the real ones follow, has no serial/parallel branch, and produces
no `log.*` files.

---

## 1. Which definition executes depends on how the case is addressed

`core/runtime/registry.py::resolve_entry`:

- name matches the registry -> `source_type: spec_factory`; the factory's
  `workflow_dag` runs. The resolution dict carries **no** `workflow_dag` key,
  so `_with_entry_metadata` leaves the factory DAG intact. An on-disk contract
  is **not consulted**.
- otherwise -> `source_type: filesystem_case`; `_with_entry_metadata`
  (registry.py:370) overwrites `metadata["workflow_dag"]` from the on-disk
  contract unconditionally — *"On-disk steps win; overwrite spec-factory
  default."*

Four directories are addressable both ways:

```
manufacturedBidomain             <-> manufacturedSolutions/bidomain
manufacturedMonodomainPseudoECG  <-> manufacturedSolutions/monodomainPseudoECG
niederer2012                     <-> NiedererEtAl2011/NiedererEtAl2011verification
singleCell                       <-> electrophysiologyProtocols/singleCell
```

Measured, same case directory, addressed both ways:

| case | by registered name | by filesystem path | agree? |
|---|---|---|---|
| `manufacturedBidomain` | blockMesh, decomposePar -force, mpirun -np 6 cardiacFoam -parallel, reconstructPar | blockMesh, cardiacFoam | **no** |
| `manufacturedMonodomainPseudoECG` | same 4-step parallel | blockMesh, cardiacFoam | **no** |
| `singleCell` | cardiacFoam | blockMesh, cardiacFoam | **no** |
| `niederer2012` | blockMesh, cardiacFoam, postProcess x2 | identical | yes |

**This is the defect.** Three of four dual-addressable cases execute
differently depending on how the entry was spelled. Nothing warns; neither
side is marked authoritative.

`mpirun -np N` versus `runParallel` is *not* part of it —
`core/runtime/parallel_execution.py` reads `numberOfSubdomains` from
`decomposeParDict`, exactly as `runParallel` does, so the rank count is
derived on both sides, not hardcoded.

---

## 2. Neither definition is reliably the correct one

This rules out the obvious fix. It is not that factories are current and
contracts stale. **They are stale in opposite directions:**

- `manufacturedBidomain`, `manufacturedMonodomainPseudoECG`: the **factory**
  is richer. Parallel execution was added there; the contract never caught up.
- `singleCell`: the **contract** is correct. The factory's DAG is a single
  unconditional `cardiacFoam` step with no mesh, but `constant/polyMesh` has
  **zero tracked files** — the mesh is not committed, so on a clean clone it
  does not exist. (The case's `Allrun` also runs `blockMesh` first, which
  corroborates the intent without being the standard.)

So "declare one side canonical" resolves nothing. Each divergent case needs
adjudicating on its own merits.

*Not verified:* whether `cardiacFoam` actually fails without a mesh for the 0D
solver. That needs a sourced OpenFOAM. The discrepancy is verified; the
runtime consequence is not.

---

## 3. Structural gaps behind the divergence

### 3.1 The DAG has no "once per study" scope

Every step in a `workflow_dag` runs per case. For a sweep over ionic models on
a fixed one-cell mesh, `blockMesh` re-runs for every case though the mesh is
identical each time. There is no way to say "this step belongs to the case
family, not the case."

This is why mesh steps drift in and out of definitions: omitting one is free
on a developer's disk where the mesh already exists, and fatal on a clean
clone. `singleCell` is exactly that failure.

### 3.2 One filename, two lifecycles

`case_root/workflow_contract.json` is an authored, committed input.
`plugins/cardiacfoam/sweep.py:125` **generates** a file of the same name into
each materialized sweep case. Same name, same reader, opposite lifecycles — a
materialized case is indistinguishable from an authored one by inspection.

This is the `run_manifest.json` confusion (retired in 5cfb8f46) still live in
a different file.

### 3.3 `_ENTRY_HINTS` is dead

`registry.py:68` declares `_ENTRY_HINTS: dict[str, dict[str, object]] = {}`
and nothing populates it. It is read at lines 197 and 199 to guard a branch
setting `entry_kind`, `source_type`, `workflow_family` and `is_runnable` — a
branch that is unreachable. It reads like a live override mechanism. Same
species as `TutorialSpec.run_case`, removed in 63a1c6e6.

---

## 4. What "one source of truth" has to decide

Stated as questions, because each has a real trade-off.

1. **Do registered tutorials keep contracts at all?** The four dual-addressable
   cases are the entire problem surface. If a registered tutorial's contract
   were derived from its factory — or simply absent — divergence becomes
   impossible rather than merely discouraged. The cost: a case folder can then
   only be run by path if it is *not* registered, which changes what
   `filesystem_case` means.
2. **Must name-addressing and path-addressing agree?** If yes, resolution can
   no longer depend on how the entry was spelled. This is the invariant worth
   asserting in a test regardless of which fix is chosen.
3. **Do steps get a scope?** Per-study versus per-case (§3.1). Without it,
   mesh steps keep drifting.
4. **What happens to the three divergent cases?** Each needs a decision about
   which behaviour was intended — physics review, not refactoring. It should
   not ride along with a mechanical change.
5. **Does the sweep materializer keep authoring `Allrun`?** (§0.2) If
   driverFOAM does not own `Allrun`, the generated one is either a layering
   violation to remove, or a deliberate exception to document.

---

## 5. Relationship to the sweep work

The isolated-workspace design (architecture-pass spec §1, option B) assumes a
case has one definition of how to run. For three of the four dual-addressable
cases, it has two that disagree. Building isolated workspaces first would mean
building them around an unresolved ambiguity about what executes inside them.

Recommend settling §4.1 and §4.2 before option B starts. §4.4 can proceed
independently and in parallel, since it is case-by-case physics review.
