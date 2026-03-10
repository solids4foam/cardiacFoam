# ECG ElectroModels Refactor — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce `ecgElectroBase` between `monoDomainElectro` and the ECG subclasses, split `greensFunctionECGElectro` into `pseudoECGElectro` and `BEMECGElectro`, eliminate 4 copies of the conductivity tensor init pattern, and add a `postProcess` flag that reads Vm from disk instead of solving the monodomain PDE.

**Architecture:** `ecgElectroBase` extends `monoDomainElectro` using the existing protected `(const word& type, Time& runTime, const word& region)` constructor chain. Shared ECG state (`Gi_`, `sigmaT_`, `postProcess_`) lives in the base; subclasses implement only `computeECG()`. `monoDomainElectro` gains two protected members: `initialiseConductivityTensor()` (static utility, replaces 4 duplicate patterns — used for both `conductivity_` and `Gi_`) and `readVm()` (called by the base in postProcess mode). The `Gi` tensor stays as a full `volTensorField` initialized from the dict (no scalar simplification — deferred).

**Tech Stack:** OpenFOAM C++ (wmake build system), existing electroModel runtime-selection pattern (`defineTypeNameAndDebug`, `addToRunTimeSelectionTable`), `volTensorField`, `dimensionedScalar`, `IOobject::READ_IF_PRESENT`.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.H` | Declare `initialiseConductivityTensor()` (protected static) and `readVm()` (protected) |
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.C` | Implement both; refactor `initialiseConductivity()` to delegate to the utility |
| **Create** | `src/electroModels/ecgElectroBase/ecgElectroBase.H` | Abstract base: `Gi_` (volTensorField), `sigmaT_`, `postProcess_`; pure virtual `computeECG()`; `evolve()` mode switch |
| **Create** | `src/electroModels/ecgElectroBase/ecgElectroBase.C` | Constructor (init `Gi_` via `initialiseConductivityTensor`), `evolve()`, `read()` |
| **Create** | `src/electroModels/pseudoECGElectro/pseudoECGElectro.H` | Gima-Rudy dipole subclass |
| **Create** | `src/electroModels/pseudoECGElectro/pseudoECGElectro.C` | `computeECG()`: dipole integral, electrode output |
| **Create** | `src/electroModels/BEMECGElectro/BEMECGElectro.H` | Green's function subclass (replaces `greensFunctionECGElectro`) |
| **Create** | `src/electroModels/BEMECGElectro/BEMECGElectro.C` | `computeECG()`: Is_ computation, Green's function integral only, torso VTK — no pseudo-ECG output |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.H` | Extend `ecgElectroBase`; remove `Gi_`, `sigmaT_`; add `sigmaE_`, `Ge_` |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.C` | Remove duplicate init; use `Gi()` from base; add `computeECG()` |
| Delete | `src/electroModels/greensFunctionECGElectro/` | Replaced by `BEMECGElectro` + `pseudoECGElectro` |
| Modify | `src/electroModels/Make/files` | Add new classes; remove `greensFunctionECGElectro` |
| Modify | `tutorials/ECG/constant/electroProperties` | `greensFunctionECGElectro` → `BEMECGElectro`; coeff dict renamed only — `Gi` tensor unchanged |

---

## Chunk 1: `monoDomainElectro` Infrastructure

**Goal:** Extract the conductivity tensor init pattern into a reusable static utility; add `readVm()` for postProcess mode. Existing behaviour is unchanged — `NiedererEtAl2012` must still pass.

### Task 1: Add `initialiseConductivityTensor` to `monoDomainElectro`

**Files:**
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.H`
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.C`

- [ ] **Step 1: Add declaration to `.H`**

  In the `protected` section of `monoDomainElectro`, after the existing `conductivity()` accessor, add:

  ```cpp
  //- Initialise a conductivity tensor field: try disk first, fall back to dict.
  //  Replaces the repeated initialiseConductivity() pattern in ECG subclasses.
  static tmp<volTensorField> initialiseConductivityTensor
  (
      const word& fieldName,
      const dictionary& dict,
      const fvMesh& mesh
  );
  ```

- [ ] **Step 2: Implement in `.C`**

  Add before `monoDomainElectro::monoDomainElectro(...)`:

  ```cpp
  Foam::tmp<Foam::volTensorField>
  Foam::monoDomainElectro::initialiseConductivityTensor
  (
      const word& fieldName,
      const dictionary& dict,
      const fvMesh& mesh
  )
  {
      IOobject io
      (
          fieldName,
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      );

      auto tresult = tmp<volTensorField>::New
      (
          IOobject
          (
              fieldName,
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedTensor(fieldName, dimConductivity, tensor::zero)
      );

      if (io.headerOk())
      {
          tresult.ref() = volTensorField(io, mesh);
          Info<< "Read " << fieldName << " from disk" << nl;
      }
      else
      {
          const dimensionedSymmTensor symmSigma(fieldName, dimConductivity, dict);
          tresult.ref() = dimensionedTensor
          (
              fieldName, dimConductivity, symmSigma.value() & tensor::I
          );
          Info<< fieldName << " = " << symmSigma << nl;
      }

      return tresult;
  }
  ```

  Note: `dimConductivity` may not exist — check the existing `initialiseConductivity()` for the correct dimensions used there (e.g. `dimless/dimLength/dimTime` or explicit `dimensionSet(-1,-3,3,0,0,2,0)`). Use exactly what the existing code uses.

- [ ] **Step 3: Refactor `initialiseConductivity()` to delegate**

  In `monoDomainElectro.C`, replace the body of `initialiseConductivity()` with:

  ```cpp
  tmp<volTensorField> monoDomainElectro::initialiseConductivity() const
  {
      return initialiseConductivityTensor
      (
          "conductivity",
          cardiacProperties_,
          mesh()
      );
  }
  ```

  The function signature and return type stay unchanged. Only the body changes.

- [ ] **Step 4: Build and verify no regressions**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles with no errors or warnings related to the change.

  Then verify the NiedererEtAl2012 tutorial still runs:

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials/NiedererEtAl2012
  ./Allrun
  ```

  Expected: completes without error; `postProcessing/` directory populated as before.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/monoDomainElectro/monoDomainElectro.H \
          src/electroModels/monoDomainElectro/monoDomainElectro.C
  git commit -m "refactor(monoDomainElectro): extract initialiseConductivityTensor static utility"
  ```

---

### Task 2: Add `readVm()` to `monoDomainElectro`

**Files:**
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.H`
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.C`

- [ ] **Step 1: Add declaration to `.H`**

  In the `protected` section, near `initialiseConductivityTensor`:

  ```cpp
  //- Re-read Vm from the current time directory (used by ecgElectroBase
  //  in postProcess mode to skip the monodomain solve).
  void readVm();
  ```

- [ ] **Step 2: Implement in `.C`**

  ```cpp
  void Foam::monoDomainElectro::readVm()
  {
      Vm_.read();
  }
  ```

- [ ] **Step 3: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly.

- [ ] **Step 4: Commit**

  ```bash
  git add src/electroModels/monoDomainElectro/monoDomainElectro.H \
          src/electroModels/monoDomainElectro/monoDomainElectro.C
  git commit -m "feat(monoDomainElectro): add readVm() for postProcess mode"
  ```

---

## Chunk 2: `ecgElectroBase`

**Goal:** New abstract class that owns shared ECG state and the evolve() mode switch. Does not register a runtime selector (it's abstract — only concrete subclasses do).

### Task 3: Create `ecgElectroBase`

**Files:**
- Create: `src/electroModels/ecgElectroBase/ecgElectroBase.H`
- Create: `src/electroModels/ecgElectroBase/ecgElectroBase.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `ecgElectroBase.H`**

  ```cpp
  #ifndef ecgElectroBase_H
  #define ecgElectroBase_H

  #include "monoDomainElectro.H"

  namespace Foam
  {

  /*---------------------------------------------------------------------------*\
                        Class ecgElectroBase Declaration
  \*---------------------------------------------------------------------------*/

  class ecgElectroBase
  :
      public monoDomainElectro
  {
      // Private data

          //- Intracellular conductivity tensor (time-invariant, set once in ctor).
          //  Initialized via initialiseConductivityTensor("Gi", ecgDict, mesh()).
          //  Gi tensor simplification to scalar (sigmaI * conductivity) is deferred.
          volTensorField Gi_;

          //- Isotropic torso / infinite-medium conductivity (constant)
          const dimensionedScalar sigmaT_;

          //- Skip monodomain solve and read Vm from disk instead
          const Switch postProcess_;


      // Private member functions

          //- ECG-type-specific computation; called every timestep after Vm is ready
          virtual void computeECG() = 0;


  public:

      TypeName("ecgElectroBase");


      // Constructors

          ecgElectroBase
          (
              const word& type,
              Time& runTime,
              const word& region = dynamicFvMesh::defaultRegion
          );


      // Destructor
      virtual ~ecgElectroBase() = default;


      // Member functions

          // Access

              const volTensorField& Gi()     const { return Gi_; }
              const dimensionedScalar& sigmaT() const { return sigmaT_; }
              bool postProcess()             const { return postProcess_; }

          // Evolution

              //- Advance one timestep: mode-switch then call computeECG()
              virtual bool evolve();

          // I/O

              virtual bool read();
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace Foam

  #endif
  ```

- [ ] **Step 2: Write `ecgElectroBase.C`**

  ```cpp
  #include "ecgElectroBase.H"

  // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  namespace Foam
  {
      defineTypeNameAndDebug(ecgElectroBase, 0);
      // No addToRunTimeSelectionTable — ecgElectroBase is abstract.
      // Only concrete subclasses register with the selector.
  }


  // * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

  Foam::ecgElectroBase::ecgElectroBase
  (
      const word& type,
      Time& runTime,
      const word& region
  )
  :
      monoDomainElectro(type, runTime, region),
      Gi_
      (
          initialiseConductivityTensor
          (
              "Gi",
              electroProperties().subDict(type + "Coeffs").subDict("ECG"),
              mesh()
          )
      ),
      sigmaT_
      (
          "sigmaT",
          electroProperties().subDict(type + "Coeffs").subDict("ECG")
      ),
      postProcess_
      (
          electroProperties().subDict(type + "Coeffs")
              .lookupOrDefault<Switch>("postProcess", false)
      )
  {
      Info<< "ECGElectroBase:" << nl
          << "  Gi      read from ECG subdict or disk" << nl
          << "  sigmaT  = " << sigmaT_.value() << " S/m" << nl
          << "  postProcess = " << postProcess_ << nl;
  }


  // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

  bool Foam::ecgElectroBase::evolve()
  {
      bool ok = true;

      if (!postProcess_)
      {
          ok = monoDomainElectro::evolve();
      }
      else
      {
          readVm();
      }

      computeECG();

      return ok;
  }


  bool Foam::ecgElectroBase::read()
  {
      // sigmaI_ and sigmaT_ are time-invariant constants — no re-read required.
      return monoDomainElectro::read();
  }
  ```

- [ ] **Step 3: Add to `Make/files`**

  In `src/electroModels/Make/files`, add after `monoDomainElectro/monoDomainElectro.C`:

  ```
  ecgElectroBase/ecgElectroBase.C
  ```

- [ ] **Step 4: Build**

  Always run `wmake` from the library root so `wmakeLnInclude` runs and `lnInclude/` symlinks are kept current. New class headers must be discoverable by sibling directories.

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly. No tests yet — the class is abstract.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/ecgElectroBase/ \
          src/electroModels/Make/files
  git commit -m "feat: add ecgElectroBase abstract class with postProcess mode switch"
  ```

---

## Chunk 3: `pseudoECGElectro` and `BEMECGElectro`

**Goal:** Two concrete ECG subclasses. `pseudoECGElectro` holds the Gima-Rudy dipole logic extracted from `greensFunctionECGElectro`. `BEMECGElectro` holds the Green's function / current-source logic. Tutorial ECG switched to `BEMECGElectro`.

### Task 4: Create `pseudoECGElectro`

**Files:**
- Create: `src/electroModels/pseudoECGElectro/pseudoECGElectro.H`
- Create: `src/electroModels/pseudoECGElectro/pseudoECGElectro.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `pseudoECGElectro.H`**

  ```cpp
  #ifndef pseudoECGElectro_H
  #define pseudoECGElectro_H

  #include "ecgElectroBase.H"
  #include "OFstream.H"

  namespace Foam
  {

  class pseudoECGElectro
  :
      public ecgElectroBase
  {
      // Private data

          wordList          electrodeNames_;
          List<vector>      electrodePositions_;
          autoPtr<OFstream> outputPtr_;    // pseudoECG.dat


      // Private member functions

          void readElectrodes(const dictionary& ecgDict);
          void writeHeader();
          virtual void computeECG();


  public:

      TypeName("pseudoECGElectro");


      // Constructors

          pseudoECGElectro
          (
              Time& runTime,
              const word& region = dynamicFvMesh::defaultRegion
          );

      virtual ~pseudoECGElectro() = default;
  };

  } // End namespace Foam

  #endif
  ```

- [ ] **Step 2: Write `pseudoECGElectro.C`**

  Copy the electrode-reading and pseudo-ECG integral from `greensFunctionECGElectro.C`. Key points:

  - `readElectrodes()`: identical to existing `greensFunctionECGElectro::readElectrodes()` — reads `ECG.electrodes` subdict.
  - `writeHeader()`: writes column headers to `pseudoECG.dat`.
  - `computeECG()`:

  ```cpp
  void Foam::pseudoECGElectro::computeECG()
  {
      const auto tgradVm = fvc::grad(Vm());
      const auto& gradVm = tgradVm();

      const scalarField& Vols = mesh().V();
      const vectorField& C    = mesh().C();
      const label nPoints     = electrodePositions_.size();

      scalarList localSums(nPoints, 0.0);

      forAll(C, cI)
      {
          forAll(electrodePositions_, pI)
          {
              const vector r_vec = C[cI] - electrodePositions_[pI];
              const scalar r     = Foam::mag(r_vec) + SMALL;
              const vector dipole = (Gi()[cI] & gradVm[cI]) * Vols[cI];
              localSums[pI] -= (dipole & r_vec) / (r*r*r);
          }
      }

      Foam::reduce(localSums, sumOp<scalarList>());

      if (runTime().outputTime() && Pstream::master())
      {
          *outputPtr_ << runTime().timeOutputValue();
          for (const scalar v : localSums)
              *outputPtr_ << tab << v;
          *outputPtr_ << nl;
      }
  }
  ```

  Use `runTime()` (from `electroModel`) rather than `mesh().time()` — they resolve to the same object but `runTime()` matches the style of the rest of the codebase.

  The factor `1/(4π sigmaT)` was absent in the existing pseudo-ECG output (it was not applied there — check the existing code and match exactly). Look at `greensFunctionECGElectro.C` lines ~360–390 for the exact formula.

  Constructor: reads `ECG` subdict, calls `readElectrodes()`, opens `postProcessing/pseudoECG.dat`, calls `writeHeader()`.

- [ ] **Step 3: Add to `Make/files`**

  ```
  pseudoECGElectro/pseudoECGElectro.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/pseudoECGElectro/ \
          src/electroModels/Make/files
  git commit -m "feat: add pseudoECGElectro (Gima-Rudy dipole) extending ecgElectroBase"
  ```

---

### Task 5: Create `BEMECGElectro`

**Files:**
- Create: `src/electroModels/BEMECGElectro/BEMECGElectro.H`
- Create: `src/electroModels/BEMECGElectro/BEMECGElectro.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `BEMECGElectro.H`**

  ```cpp
  #ifndef BEMECGElectro_H
  #define BEMECGElectro_H

  #include "ecgElectroBase.H"
  #include "OFstream.H"
  #include "triSurface.H"

  namespace Foam
  {

  class BEMECGElectro
  :
      public ecgElectroBase
  {
      // Private data

          //- Current source density: Is = -div(Gi & grad(Vm))
          volScalarField Is_;

          //- Electrode names and positions
          wordList     electrodeNames_;
          List<vector> electrodePositions_;

          //- Torso surface (optional)
          autoPtr<triSurface> torsoSurfacePtr_;
          pointField          torsoFaceCentres_;
          fileName            torsoVtkDir_;

          //- Output stream
          autoPtr<OFstream> outputPtr_;   // ECG.dat only — pseudo-ECG belongs to pseudoECGElectro


      // Private member functions

          void readElectrodes(const dictionary& ecgDict);
          void loadTorsoSurface(const dictionary& ecgDict);
          void writeHeader();
          void writeTorsoVtk(const scalarList& phiGreens) const;
          virtual void computeECG();


  public:

      TypeName("BEMECGElectro");


      // Constructors

          BEMECGElectro
          (
              Time& runTime,
              const word& region = dynamicFvMesh::defaultRegion
          );

      virtual ~BEMECGElectro() = default;
  };

  } // End namespace Foam

  #endif
  ```

- [ ] **Step 2: Write `BEMECGElectro.C`**

  This is a direct refactor of `greensFunctionECGElectro.C`. Key differences from the original:

  1. Base class: `ecgElectroBase` instead of `monoDomainElectro`
  2. **Remove** `Gi_` member and `initialiseGi()` — use `Gi()` from base
  3. **Remove** `sigmaT_` member — use `sigmaT()` from base
  4. **Remove** `outputPseudoPtr_` and `writePseudoHeader()` — pseudo-ECG output is `pseudoECGElectro`'s responsibility
  5. `computeECG()` replaces the ECG portion of `evolve()` — the monodomain call is gone

  Constructor: copy from `greensFunctionECGElectro`, remove `Gi_`/`sigmaT_` init and the pseudo-ECG output stream setup.

  `computeECG()`: copy the ECG block from `greensFunctionECGElectro::evolve()` (lines 302–409 of original):
  - Keep the Green's function integral (`localSumsGreens`) and `ECG.dat` output
  - **Drop** the pseudo-ECG integral (`localSumsPseudo`) and `pseudoECG.dat` output entirely
  - Replace `Gi_` with `Gi()` and `sigmaT_` with `sigmaT()`

  `writeTorsoVtk` takes only `phiGreens` (no `phiPseudo`).

  Methods `readElectrodes`, `loadTorsoSurface`, `writeHeader` copy verbatim from `greensFunctionECGElectro.C`.

- [ ] **Step 3: Add to `Make/files`**

  ```
  BEMECGElectro/BEMECGElectro.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/BEMECGElectro/ \
          src/electroModels/Make/files
  git commit -m "feat: add BEMECGElectro (Green's function ECG) extending ecgElectroBase"
  ```

---

### Task 6: Update tutorial and verify `BEMECGElectro` end-to-end

**Files:**
- Modify: `tutorials/ECG/constant/electroProperties`

- [ ] **Step 1: Update `electroProperties`**

  Two changes only — everything else stays verbatim:
  1. `electroModel` value: `greensFunctionECGElectro` → `BEMECGElectro`
  2. Coeff dict name: `greensFunctionECGElectroCoeffs` → `BEMECGElectroCoeffs`

  The `ECG {}` subdict content is **unchanged** — `Gi` stays as the full anisotropic tensor:

  ```
  electroModel  BEMECGElectro;

  BEMECGElectroCoeffs
  {
      chi   [0 -1 0 0 0 0 0] 140000;
      Cm    [-1 -4 4 0 0 2 0] 0.01;
      ionicModel    BuenoOrovio;
      tissue        epicardialCells;
      solutionAlgorithm    explicit;
      // ... all other monodomain keys unchanged ...

      ECG
      {
          Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);  // unchanged
          sigmaT  [-1 -3 3 0 0 2 0] 0.24725;                      // unchanged
          electrodes
          {
              V1 (-26.1184 -283.543 -70.0558);
              // ... all electrodes unchanged ...
          }
      }
  }
  ```

- [ ] **Step 2: Run ECG tutorial**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials/ECG
  ./Allclean && ./Allrun
  ```

  Expected:
  - Solver runs without error
  - `postProcessing/ECG.dat` and `postProcessing/pseudoECG.dat` populated with time series
  - ECG signals comparable to previous results

- [ ] **Step 3: Commit**

  ```bash
  git add tutorials/ECG/constant/electroProperties
  git commit -m "feat(tutorials/ECG): switch to BEMECGElectro with sigmaI scalar"
  ```

---

## Chunk 4: `pdeECGElectro` Refactor and Cleanup

**Goal:** Migrate `pdeECGElectro` to extend `ecgElectroBase`, eliminate its duplicate members, then delete `greensFunctionECGElectro`.

### Task 7: Refactor `pdeECGElectro`

**Files:**
- Modify: `src/electroModels/pdeECGElectro/pdeECGElectro.H`
- Modify: `src/electroModels/pdeECGElectro/pdeECGElectro.C`

- [ ] **Step 1: Update `.H`**

  Change base class from `monoDomainElectro` to `ecgElectroBase`:
  ```cpp
  #include "ecgElectroBase.H"

  class pdeECGElectro : public ecgElectroBase
  ```

  **Remove** private members:
  - `volTensorField Gi_` — now in `ecgElectroBase`
  - `dimensionedScalar sigmaT_` — now in `ecgElectroBase`

  **Add** private members:
  ```cpp
  const dimensionedScalar sigmaE_;   // extracellular scalar: Ge = sigmaE_ * conductivity()
  volTensorField          Ge_;
  ```

  Keep unchanged: `phiE_`, `GiPlusGe_`, `torsoMeshPtr_`, `phiTPtr_`, `heartPatchName_`, `torsoPatchName_`.

  **Remove** private methods:
  - `tmp<volTensorField> initialiseGi() const`

  **Add** pure virtual implementation:
  ```cpp
  virtual void computeECG();
  ```

  **Remove** override of `evolve()` — base class handles it now.

- [ ] **Step 2: Update `.C` constructor**

  Change base call from `monoDomainElectro(type, runTime, region)` to `ecgElectroBase(type, runTime, region)`.

  **Remove**:
  - `Gi_` initialiser (use `Gi()` from base)
  - `sigmaT_` initialiser (use `sigmaT()` from base)
  - `initialiseGi()` method body

  **Add** in initialiser list:
  ```cpp
  sigmaE_
  (
      "sigmaE",
      electroProperties().subDict(type + "Coeffs").subDict("ECG")
  ),
  Ge_
  (
      IOobject("Ge", runTime.timeName(), mesh(), IOobject::NO_READ, IOobject::NO_WRITE),
      sigmaE_ * conductivity()
  ),
  ```

  Update `GiPlusGe_` initialisation to use `Gi()` from base:
  ```cpp
  GiPlusGe_(IOobject(...), Gi() + Ge_),
  ```

- [ ] **Step 3: Rename `evolve()` body to `computeECG()`**

  The current `pdeECGElectro::evolve()`:
  1. calls `monoDomainElectro::evolve()` — **remove this line** (base handles it)
  2. solves Poisson on heart — **keep**
  3. updates torso BCs — **keep**
  4. solves Laplace on torso — **keep**
  5. writes at `outputTime()` — **keep**

  Replace all references to `Gi_` with `Gi()`, `sigmaT_` with `sigmaT()`.

  Rename method to `computeECG()`.

- [ ] **Step 4: Update `electroProperties` for pdeECGElectro (if a test case exists)**

  If a tutorial uses `pdeECGElectro`, update its dict to use `sigmaI` and `sigmaE` scalars instead of `Gi`/`Ge` tensors. Check for any test case under `tutorials/` first:

  ```bash
  grep -r pdeECGElectro /Users/simaocastro/cardiacFoamv2/tutorials/
  ```

- [ ] **Step 5: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly.

- [ ] **Step 6: Commit**

  ```bash
  git add src/electroModels/pdeECGElectro/
  git commit -m "refactor(pdeECGElectro): extend ecgElectroBase, remove duplicate Gi/sigmaT"
  ```

---

### Task 8: Delete `greensFunctionECGElectro` and final cleanup

**Files:**
- Delete: `src/electroModels/greensFunctionECGElectro/`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Remove `greensFunctionECGElectro` from `Make/files`**

  Delete the line:
  ```
  greensFunctionECGElectro/greensFunctionECGElectro.C
  ```

- [ ] **Step 2: Delete the directory**

  ```bash
  rm -rf /Users/simaocastro/cardiacFoamv2/src/electroModels/greensFunctionECGElectro
  ```

- [ ] **Step 3: Build clean**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly. No references to `greensFunctionECGElectro` should remain.

- [ ] **Step 4: Full regression**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials
  ./Alltest-regression
  ```

  Expected: all tutorials pass.

- [ ] **Step 5: Commit**

  ```bash
  git add -A src/electroModels/
  git commit -m "chore: remove greensFunctionECGElectro (replaced by BEMECGElectro + pseudoECGElectro)"
  ```

---

## Deferred / Known Limitations

1. **`Gi` tensor → scalar simplification**: `Gi_` stays as a full anisotropic `volTensorField` initialized from the dict. The physical relationship `Gi = sigmaI * conductivity()` is noted in the design spec but not implemented here — a separate task once the anisotropy handling is better understood.

2. **`postProcess` time-looping**: `readVm()` re-reads Vm from `mesh().time().timeName()`. The solver binary must advance `runTime` through existing time directories for this to be useful. Integration with the solver binary or a dedicated utility is **deferred** — the flag is wired up but standalone postProcess invocation is not tested in this plan.
