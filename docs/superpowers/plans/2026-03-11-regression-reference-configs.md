# Regression Reference Configs Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ensure regression tests always run with the canonical `electroProperties` config, regardless of user edits, by storing a `.reference` copy and restoring it before each test run.

**Architecture:** Each tutorial keeps `constant/electroProperties.reference` as a locked copy of the config the reference values were captured with. `runRegressionTest.sh` copies all `constant/*.reference` files to their canonical names before invoking the solver.

**Tech Stack:** Bash, OpenFOAM dict files.

---

## Chunk 1: singleCell

**Files:**
- Create: `tutorials/singleCell/constant/electroProperties.reference`
- Modify: `tutorials/singleCell/runRegressionTest.sh`

### Task 1: Create the reference config for singleCell

- [ ] **Step 1: Copy the canonical electroProperties**

```bash
cp tutorials/singleCell/constant/electroProperties \
   tutorials/singleCell/constant/electroProperties.reference
```

- [ ] **Step 2: Verify the file exists and matches**

```bash
diff tutorials/singleCell/constant/electroProperties \
     tutorials/singleCell/constant/electroProperties.reference
```
Expected: no output (identical files).

- [ ] **Step 3: Commit**

```bash
git add tutorials/singleCell/constant/electroProperties.reference
git commit -m "test(singleCell): add reference electroProperties for regression test"
```

---

### Task 2: Add restoreReferenceConfigs to singleCell runRegressionTest.sh

- [ ] **Step 1: Read the current script**

Open `tutorials/singleCell/runRegressionTest.sh` and locate the line `. $WM_PROJECT_DIR/bin/tools/RunFunctions`.

- [ ] **Step 2: Add the restore function after the source line**

Insert the following block immediately after the `. $WM_PROJECT_DIR/bin/tools/RunFunctions` line:

```bash
restoreReferenceConfigs()
{
    for refFile in constant/*.reference
    do
        [ -f "${refFile}" ] || continue
        target="${refFile%.reference}"
        echo "Restoring ${target} from ${refFile}"
        cp "${refFile}" "${target}"
    done
}
```

- [ ] **Step 3: Call restoreReferenceConfigs before the solver**

Add a call to `restoreReferenceConfigs` immediately before the `runApplication cardiacFoam` line:

```bash
restoreReferenceConfigs

runApplication cardiacFoam
```

- [ ] **Step 4: Verify the script looks correct**

The relevant section of `runRegressionTest.sh` should read:

```bash
. $WM_PROJECT_DIR/bin/tools/RunFunctions

restoreReferenceConfigs()
{
    for refFile in constant/*.reference
    do
        [ -f "${refFile}" ] || continue
        target="${refFile%.reference}"
        echo "Restoring ${target} from ${refFile}"
        cp "${refFile}" "${target}"
    done
}

checkSingleCellResults()
{
    ...
}

restoreReferenceConfigs

runApplication cardiacFoam
```

- [ ] **Step 5: Smoke-test the restore logic in isolation**

Modify `constant/electroProperties` by changing `ionicModel` to something invalid, then confirm the restore function overwrites it:

```bash
cd tutorials/singleCell
echo "ionicModel broken;" > constant/electroProperties   # break it
bash -c '. $WM_PROJECT_DIR/bin/tools/RunFunctions
for refFile in constant/*.reference; do
    [ -f "${refFile}" ] || continue
    target="${refFile%.reference}"
    cp "${refFile}" "${target}"
done'
grep ionicModel constant/electroProperties
```

Expected output: `ionicModel    BuenoOrovio;` (restored from `.reference`).

- [ ] **Step 6: Commit**

```bash
git add tutorials/singleCell/runRegressionTest.sh
git commit -m "test(singleCell): restore reference configs before regression run"
```

---

## Chunk 2: NiedererEtAl2012

**Files:**
- Create: `tutorials/NiedererEtAl2012/constant/electroProperties.reference`
- Modify: `tutorials/NiedererEtAl2012/runRegressionTest.sh`

### Task 3: Create the reference config for NiedererEtAl2012

- [ ] **Step 1: Copy the canonical electroProperties**

```bash
cp tutorials/NiedererEtAl2012/constant/electroProperties \
   tutorials/NiedererEtAl2012/constant/electroProperties.reference
```

- [ ] **Step 2: Verify the file exists and matches**

```bash
diff tutorials/NiedererEtAl2012/constant/electroProperties \
     tutorials/NiedererEtAl2012/constant/electroProperties.reference
```
Expected: no output (identical files).

- [ ] **Step 3: Commit**

```bash
git add tutorials/NiedererEtAl2012/constant/electroProperties.reference
git commit -m "test(NiedererEtAl2012): add reference electroProperties for regression test"
```

---

### Task 4: Add restoreReferenceConfigs to NiedererEtAl2012 runRegressionTest.sh

- [ ] **Step 1: Read the current script**

Open `tutorials/NiedererEtAl2012/runRegressionTest.sh` and locate the line `. $WM_PROJECT_DIR/bin/tools/RunFunctions`.

- [ ] **Step 2: Add the restore function after the source line**

Insert the following block immediately after the `. $WM_PROJECT_DIR/bin/tools/RunFunctions` line:

```bash
restoreReferenceConfigs()
{
    for refFile in constant/*.reference
    do
        [ -f "${refFile}" ] || continue
        target="${refFile%.reference}"
        echo "Restoring ${target} from ${refFile}"
        cp "${refFile}" "${target}"
    done
}
```

- [ ] **Step 3: Call restoreReferenceConfigs before the solver**

Add a call to `restoreReferenceConfigs` immediately before the `runApplication blockMesh` line (blockMesh runs before the solver, but the configs must be in place before anything runs):

```bash
restoreReferenceConfigs

runApplication blockMesh
```

- [ ] **Step 4: Verify the script looks correct**

The top of the script (after the source line) should read:

```bash
. $WM_PROJECT_DIR/bin/tools/RunFunctions

restoreReferenceConfigs()
{
    for refFile in constant/*.reference
    do
        [ -f "${refFile}" ] || continue
        target="${refFile%.reference}"
        echo "Restoring ${target} from ${refFile}"
        cp "${refFile}" "${target}"
    done
}

checkSmokeCheckResults()
{
    ...
}

restoreReferenceConfigs

runApplication blockMesh
```

- [ ] **Step 5: Smoke-test the restore logic in isolation**

```bash
cd tutorials/NiedererEtAl2012
echo "electroModel broken;" > constant/electroProperties   # break it
bash -c '. $WM_PROJECT_DIR/bin/tools/RunFunctions
for refFile in constant/*.reference; do
    [ -f "${refFile}" ] || continue
    target="${refFile%.reference}"
    cp "${refFile}" "${target}"
done'
grep electroModel constant/electroProperties
```

Expected output: `electroModel monoDomainElectro;` (restored from `.reference`).

- [ ] **Step 6: Commit**

```bash
git add tutorials/NiedererEtAl2012/runRegressionTest.sh
git commit -m "test(NiedererEtAl2012): restore reference configs before regression run"
```

---

## Chunk 3: Documentation note

### Task 5: Update MEMORY.md with the new convention

- [ ] **Step 1: Add a note to memory**

Add the following to `~/.claude/projects/.../memory/MEMORY.md` under a "Regression Tests" section:

```
## Regression Tests
- Reference configs stored as `constant/*.reference` in each tutorial
- `runRegressionTest.sh` calls `restoreReferenceConfigs()` before the solver — copies `.reference` → canonical name
- To recapture baselines: update `electroProperties`, run solver, update `system/*.reference` values file, then copy `electroProperties` → `electroProperties.reference`
```

- [ ] **Step 2: Commit**

```bash
git add docs/superpowers/specs/2026-03-11-regression-reference-configs-design.md \
        docs/superpowers/plans/2026-03-11-regression-reference-configs.md
git commit -m "docs: add regression reference configs spec and plan"
```
