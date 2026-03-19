# Pseudo-ECG in greensFunctionECGElectro — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `greensFunctionECGElectro` to simultaneously compute and output pseudo-ECG (dipole) potentials alongside the existing Green's function (monopole) potentials, in a single cell loop pass.

**Architecture:** Extend the existing class in-place — no new class. Both integrals share the same double loop over heart cells × evaluation points. A second output stream (`pseudoECG.dat`) is added for electrode output; `writeTorsoVtk()` gains a second scalar field (`phiPseudo`) alongside the existing `phiGreens`.

**Tech Stack:** OpenFOAM C++, wmake, `fvc::grad`, `autoPtr<OFstream>`, ASCII VTK POLYDATA.

**Spec:** [docs/superpowers/specs/2026-03-10-pseudo-ecg-design.md](../specs/2026-03-10-pseudo-ecg-design.md)

---

## Chunk 1: Header + constructor + writeTorsoVtk

### Task 1: Update the header file

**Files:**

- Modify: `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.H`

- [ ] **Step 1: Read the current header**

  Open [greensFunctionECGElectro.H](../../src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.H) and locate the private data section and the `writeTorsoVtk` declaration.

- [ ] **Step 2: Add `outputPseudoPtr_` after `outputPtr_`**

  In the private data section, after:

  ```cpp
  //- Output stream for ECG.dat (null if electrodes not configured)
  autoPtr<OFstream> outputPtr_;
  ```

  add:

  ```cpp
  //- Output stream for pseudoECG.dat (null if electrodes not configured)
  autoPtr<OFstream> outputPseudoPtr_;
  ```

- [ ] **Step 3: Update `writeTorsoVtk` declaration**

  Change:

  ```cpp
  //- Write torso surface potential as VTK ASCII POLYDATA
  void writeTorsoVtk(const scalarList& phiTorso) const;
  ```

  to:

  ```cpp
  //- Write torso surface potentials as VTK ASCII POLYDATA
  void writeTorsoVtk
  (
      const scalarList& phiGreens,
      const scalarList& phiPseudo
  ) const;
  ```

- [ ] **Step 4: Add `writePseudoHeader()` declaration in the private member functions section, after `writeHeader()`**

  ```cpp
  //- Write column header to pseudoECG.dat
  void writePseudoHeader();
  ```

- [ ] **Step 5: Verify the header compiles in isolation**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake 2>&1 | head -30
  ```

  Expected: compile error referencing `writeTorsoVtk` mismatch (`.C` still has old signature) — confirms the header change is picked up. If the build unexpectedly passes, verify that `greensFunctionECGElectro.C` still has the old `writeTorsoVtk(phiTorso)` signature and was not already modified.

---

### Task 2: Update the constructor in the `.C` file

**Files:**

- Modify: `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.C`

- [ ] **Step 1: Add `outputPseudoPtr_()` to the initialiser list**

  In the constructor initialiser list, after:

  ```cpp
  outputPtr_(),
  ```

  add:

  ```cpp
  outputPseudoPtr_(),
  ```

- [ ] **Step 2: Open `pseudoECG.dat` and write its header**

  Inside the constructor body, in the `if (hasElectrodes)` block, after:

  ```cpp
  outputPtr_.reset(new OFstream(outDir/"ECG.dat"));
  writeHeader();
  ```

  add:

  ```cpp
  outputPseudoPtr_.reset(new OFstream(outDir/"pseudoECG.dat"));
  writePseudoHeader();
  ```

- [ ] **Step 3: Implement `writePseudoHeader()` in the `.C` file**

  Add immediately after the `writeHeader()` implementation (around line 137):

  ```cpp
  void greensFunctionECGElectro::writePseudoHeader()
  {
      OFstream& os = outputPseudoPtr_.ref();
      os.setf(std::ios::scientific);
      os.precision(8);

      os << "# time";
      forAll(electrodeNames_, i)
      {
          os << "  " << electrodeNames_[i];
      }
      os << nl;
  }
  ```

---

### Task 3: Update `writeTorsoVtk()` to write two scalar fields

**Files:**

- Modify: `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.C`

- [ ] **Step 1: Update the function signature**

  Change:

  ```cpp
  void greensFunctionECGElectro::writeTorsoVtk(const scalarList& phiTorso) const
  ```

  to:

  ```cpp
  void greensFunctionECGElectro::writeTorsoVtk
  (
      const scalarList& phiGreens,
      const scalarList& phiPseudo
  ) const
  ```

- [ ] **Step 2: Replace the single CELL_DATA block with two scalar sections**

  Replace the existing `CELL_DATA` block:

  ```cpp
  os << "CELL_DATA " << nTris << "\n";
  os << "SCALARS phiT float 1\n";
  os << "LOOKUP_TABLE default\n";
  forAll(phiTorso, fI)
  {
      os << phiTorso[fI] << "\n";
  }
  ```

  with:

  ```cpp
  os << "CELL_DATA " << nTris << "\n";

  os << "SCALARS phiGreens float 1\n";
  os << "LOOKUP_TABLE default\n";
  forAll(phiGreens, fI)
  {
      os << phiGreens[fI] << "\n";
  }

  os << "SCALARS phiPseudo float 1\n";
  os << "LOOKUP_TABLE default\n";
  forAll(phiPseudo, fI)
  {
      os << phiPseudo[fI] << "\n";
  }
  ```

- [ ] **Step 3: Verify the library builds**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake 2>&1 | tail -20
  ```

  Expected: compile error in `evolve()` — old call `writeTorsoVtk(phiTorso)` now has wrong arity. This confirms the signature change propagated correctly. The `.C` file will be fixed in Chunk 2. If the build unexpectedly passes, check whether `evolve()` still contains `writeTorsoVtk(phiTorso)` — if not, someone already patched it.

- [ ] **Step 4: Commit header + constructor + writeTorsoVtk changes (Tasks 1, 2, and 3)**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2
  git add src/electroModels/greensFunctionECGElectro/
  git commit -m "feat(greensFunctionECGElectro): add outputPseudoPtr_ and two-field writeTorsoVtk"
  ```

---

## Chunk 2: Update `evolve()` and integration test

### Task 4: Rewrite `evolve()` with the combined loop

**Files:**

- Modify: `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.C`

- [ ] **Step 1: Compute `gradVm` before the cell loop**

  In `evolve()`, after:

  ```cpp
  const vectorField& Ctrs = mesh().C().primitiveField();
  ```

  add:

  ```cpp
  tmp<volVectorField> tgradVm = fvc::grad(Vm());
  const vectorField& gradVm   = tgradVm().primitiveField();
  ```

- [ ] **Step 2: Replace the single accumulator with two**

  Replace:

  ```cpp
  List<scalar> localSums(nAll, scalar(0));
  ```

  with:

  ```cpp
  List<scalar> localSumsGreens(nAll, scalar(0));
  List<scalar> localSumsPseudo(nAll, scalar(0));
  ```

- [ ] **Step 3: Replace the inner loop body**

  Replace:

  ```cpp
  forAll(Ctrs, cI)
  {
      const scalar IsV = IsI[cI]*Vols[cI];
      for (label pI = 0; pI < nAll; pI++)
      {
          const scalar r = mag(Ctrs[cI] - allPoints[pI]);
          if (r > VSMALL)
          {
              localSums[pI] += IsV/r;
          }
      }
  }
  ```

  with:

  ```cpp
  forAll(Ctrs, cI)
  {
      const scalar IsV     = IsI[cI]*Vols[cI];
      const vector gradVmV = gradVm[cI]*Vols[cI];

      for (label pI = 0; pI < nAll; pI++)
      {
          const vector r_vec = Ctrs[cI] - allPoints[pI];
          const scalar r     = mag(r_vec);
          if (r > VSMALL)
          {
              localSumsGreens[pI] += IsV / r;
              localSumsPseudo[pI] += (gradVmV & r_vec) / (r*r*r);
          }
      }
  }
  ```

- [ ] **Step 4: Update the reduce + coefficient loop**

  Replace:

  ```cpp
  for (label pI = 0; pI < nAll; pI++)
  {
      reduce(localSums[pI], sumOp<scalar>());
      localSums[pI] *= invCoeff;
  }
  ```

  with:

  ```cpp
  for (label pI = 0; pI < nAll; pI++)
  {
      reduce(localSumsGreens[pI], sumOp<scalar>());
      reduce(localSumsPseudo[pI], sumOp<scalar>());
      localSumsGreens[pI] *= invCoeff;
  }
  ```

  Note: `localSumsPseudo` is NOT multiplied by `invCoeff` — pseudo-ECG has no `sigmaT` factor.

- [ ] **Step 5: Update the electrode output block**

  Replace:

  ```cpp
  if (nE > 0)
  {
      OFstream& os = outputPtr_.ref();
      os << runTime().value();
      for (label eI = 0; eI < nE; eI++)
      {
          os << "  " << localSums[eI];
      }
      os << nl;
  }
  ```

  with:

  ```cpp
  if (nE > 0)
  {
      // Green's function → ECG.dat
      {
          OFstream& os = outputPtr_.ref();
          os << runTime().value();
          for (label eI = 0; eI < nE; eI++)
          {
              os << "  " << localSumsGreens[eI];
          }
          os << nl;
      }

      // Pseudo-ECG → pseudoECG.dat
      {
          OFstream& os = outputPseudoPtr_.ref();
          os << runTime().value();
          for (label eI = 0; eI < nE; eI++)
          {
              os << "  " << localSumsPseudo[eI];
          }
          os << nl;
      }
  }
  ```

- [ ] **Step 6: Update the torso surface output block**

  Replace:

  ```cpp
  if (nF > 0 && runTime().outputTime())
  {
      scalarList phiTorso(nF);
      for (label fI = 0; fI < nF; fI++)
      {
          phiTorso[fI] = localSums[nE + fI];
      }
      writeTorsoVtk(phiTorso);
  }
  ```

  with:

  ```cpp
  if (nF > 0 && runTime().outputTime())
  {
      scalarList phiGreens(nF);
      scalarList phiPseudo(nF);
      for (label fI = 0; fI < nF; fI++)
      {
          phiGreens[fI] = localSumsGreens[nE + fI];
          phiPseudo[fI] = localSumsPseudo[nE + fI];
      }
      writeTorsoVtk(phiGreens, phiPseudo);
  }
  ```

---

### Task 5: Build and integration test

**Files:** none new

- [ ] **Step 1: Build the library**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake 2>&1 | tail -20
  ```

  Expected: clean build, `libelectroModels.so` updated, no warnings about `localSums`.

- [ ] **Step 2: Run the ECG tutorial smoke check**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials/ECG
  ./Allrun 2>&1 | tail -30
  ```

  Expected: solver runs to completion without errors.

- [ ] **Step 3: Verify `pseudoECG.dat` is created**

  ```bash
  head -5 /Users/simaocastro/cardiacFoamv2/tutorials/ECG/postProcessing/pseudoECG.dat
  ```

  Expected: header line starting with `# time` followed by electrode names, then numeric rows.

- [ ] **Step 4: Verify `ECG.dat` still has the same format**

  ```bash
  head -5 /Users/simaocastro/cardiacFoamv2/tutorials/ECG/postProcessing/ECG.dat
  ```

  Expected: same format as before — header + numeric rows. Confirm it was not broken.

- [ ] **Step 5: Verify VTK contains both scalar fields (if torso surface configured)**

  If the ECG tutorial has a `torsoSurface` entry:

  ```bash
  grep -c "SCALARS" /Users/simaocastro/cardiacFoamv2/tutorials/ECG/postProcessing/ECG/phiT_*.vtk | head -3
  ```

  Expected: each VTK file reports `2` (one `phiGreens`, one `phiPseudo`).

- [ ] **Step 6: Commit**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2
  git add src/electroModels/greensFunctionECGElectro/
  git commit -m "feat(greensFunctionECGElectro): add pseudo-ECG dipole integral alongside Green's function"
  ```
