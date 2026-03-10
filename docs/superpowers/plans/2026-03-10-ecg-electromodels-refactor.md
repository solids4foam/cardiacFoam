# ECG ElectroModels Refactor — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce `ecgElectroBase` between `monoDomainElectro` and the ECG subclasses, split `greensFunctionECGElectro` into `pseudoECGElectro` and `BEMECGElectro`, eliminate 4 copies of the conductivity tensor init pattern, and add a `postProcess` flag that reads Vm from disk instead of solving the monodomain PDE.

**Architecture:** `ecgElectroBase` extends `monoDomainElectro` using the existing protected `(const word& type, Time& runTime, const word& region)` constructor chain. Shared ECG state (`sigmaI_`, `Gi_`, `sigmaT_`) lives in the base; subclasses implement only `computeECG()`. `monoDomainElectro` gains two protected members: `initialiseConductivityTensor()` (static utility, replaces 4 duplicate patterns) and `readVm()` (called by the base in postProcess mode).

**Tech Stack:** OpenFOAM C++ (wmake build system), existing electroModel runtime-selection pattern (`defineTypeNameAndDebug`, `addToRunTimeSelectionTable`), `volTensorField`, `dimensionedScalar`, `IOobject::READ_IF_PRESENT`.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.H` | Declare `initialiseConductivityTensor()` (protected static) and `readVm()` (protected) |
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.C` | Implement both; refactor `initialiseConductivity()` to delegate to the utility |
| **Create** | `src/electroModels/ecgElectroBase/ecgElectroBase.H` | Abstract base: `sigmaI_`, `Gi_`, `sigmaT_`, `postProcess_`; pure virtual `computeECG()`; `evolve()` mode switch |
| **Create** | `src/electroModels/ecgElectroBase/ecgElectroBase.C` | Constructor, `evolve()`, `read()` |
| **Create** | `src/electroModels/pseudoECGElectro/pseudoECGElectro.H` | Gima-Rudy dipole subclass |
| **Create** | `src/electroModels/pseudoECGElectro/pseudoECGElectro.C` | `computeECG()`: dipole integral, electrode output |
| **Create** | `src/electroModels/BEMECGElectro/BEMECGElectro.H` | Green's function subclass (replaces `greensFunctionECGElectro`) |
| **Create** | `src/electroModels/BEMECGElectro/BEMECGElectro.C` | `computeECG()`: Is_ computation, Green's function integral, torso VTK |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.H` | Extend `ecgElectroBase`; remove `Gi_`, `sigmaT_`; add `sigmaE_`, `Ge_` |
| Modify | `src/electroModels/pdeECGElectro/pdeECGElectro.C` | Remove duplicate init; use `Gi()` from base; add `computeECG()` |
| Delete | `src/electroModels/greensFunctionECGElectro/` | Replaced by `BEMECGElectro` + `pseudoECGElectro` |
| Modify | `src/electroModels/Make/files` | Add new classes; remove `greensFunctionECGElectro` |
| Modify | `tutorials/ECG/constant/electroProperties` | `greensFunctionECGElectro` → `BEMECGElectro`; `Gi` tensor → `sigmaI` scalar |

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

          //- Intracellular conductivity scalar: Gi_ = sigmaI_ * conductivity()
          const dimensionedScalar sigmaI_;

          //- Intracellular conductivity tensor (time-invariant, set once in ctor)
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

              const volTensorField& Gi() const { return Gi_; }
              const dimensionedScalar& sigmaT() const { return sigmaT_; }
              bool postProcess() const { return postProcess_; }

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
      sigmaI_
      (
          "sigmaI",
          electroProperties().subDict(type + "Coeffs").subDict("ECG")
      ),
      Gi_
      (
          IOobject
          (
              "Gi",
              runTime.timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          sigmaI_ * conductivity()
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
          << "  sigmaI  = " << sigmaI_.value() << " S/m" << nl
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

          //- Output streams
          autoPtr<OFstream> outputPtr_;       // ECG.dat
          autoPtr<OFstream> outputPseudoPtr_; // pseudoECG.dat (companion output)


      // Private member functions

          void readElectrodes(const dictionary& ecgDict);
          void loadTorsoSurface(const dictionary& ecgDict);
          void writeHeader();
          void writePseudoHeader();
          void writeTorsoVtk
          (
              const scalarList& phiGreens,
              const scalarList& phiPseudo
          ) const;
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
  4. `computeECG()` replaces the ECG portion of `evolve()` — the monodomain call is gone

  Constructor: copy from `greensFunctionECGElectro`, remove all `Gi_`/`sigmaT_` init.

  `computeECG()`: copy the ECG block from `greensFunctionECGElectro::evolve()` (lines 302–409 of original), replacing `Gi_` with `Gi()` and `sigmaT_` with `sigmaT()`.

  All other methods (`readElectrodes`, `loadTorsoSurface`, `writeHeader`, `writePseudoHeader`, `writeTorsoVtk`) copy verbatim from `greensFunctionECGElectro.C`.

  The single-pass loop combining Green's function + pseudo-ECG integrals is preserved in `computeECG()` for efficiency.

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

  Change:
  ```
  electroModel  greensFunctionECGElectro;

  greensFunctionECGElectroCoeffs
  {
      ...
      ECG
      {
          Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);
          sigmaT  [-1 -3 3 0 0 2 0] 0.24725;
          ...
      }
  }
  ```

  To:
  ```
  electroModel  BEMECGElectro;

  BEMECGElectroCoeffs
  {
      // All monodomain settings remain flat here (chi, Cm, ionicModel, etc.)
      // — copy verbatim from greensFunctionECGElectroCoeffs above
      chi   [0 -1 0 0 0 0 0] 140000;
      Cm    [-1 -4 4 0 0 2 0] 0.01;
      ionicModel    BuenoOrovio;
      tissue        epicardialCells;
      solutionAlgorithm    explicit;
      // ... (all other monodomain keys unchanged)

      ECG
      {
          sigmaI  [-1 -3 3 0 0 2 0] 0.17;     // was: Gi tensor (6 components)
          sigmaT  [-1 -3 3 0 0 2 0] 0.24725;  // unchanged
          electrodes
          {
              V1 (-26.1184 -283.543 -70.0558);
              // ... all other electrodes unchanged
          }
      }
  }
  ```

  Note: `sigmaI` is the isotropic intracellular conductivity scalar that was the diagonal of the original `Gi` tensor. Use `0.17` for the longitudinal component (the existing `Gi` was anisotropic with `(0.17 0 0  0.019 0  0.019)` — confirm with the user whether `sigmaI = 0.17` (fibre direction) or a harmonic mean should be used here).

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

## Open Questions (Confirm Before Implementing)

1. **`sigmaI` value in tutorial**: The original `Gi` was anisotropic `(0.17 0 0  0.019 0  0.019)`. The new `sigmaI` is a scalar. Confirm whether `sigmaI = 0.17` (fibre-direction value) is appropriate, or if a different representative value should be used.

2. **`postProcess` time-looping**: In postProcess mode, `readVm()` re-reads from `mesh().time().timeName()`. The solver binary must advance `runTime` through existing time directories. This integration with the solver binary or a dedicated utility is **deferred** — the flag is implemented but not the solver-level looping.

3. **`pseudoECGElectro` companion output**: Currently `greensFunctionECGElectro` outputs both `ECG.dat` (Green's) and `pseudoECG.dat` (dipole) in a single pass. After splitting, decide whether `BEMECGElectro` keeps its companion `pseudoECG.dat` output or drops it entirely (since `pseudoECGElectro` now handles it).
