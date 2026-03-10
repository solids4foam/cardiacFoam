# ECG ElectroModels Refactor — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce `ecgElectro` (concrete electroModel) + `ecgModel` (runtime-selectable, like `ionicModel`) to replace `greensFunctionECGElectro`. ECG type is selected inside the electroModel coeff dict (`ecgModel BEMECGElectro;`) and has its own subdict, exactly like `ionicModel`.

**Architecture:** `ecgElectro` extends `monoDomainElectro` and owns an `autoPtr<ecgModel>`. `ecgModel` is an abstract runtime-selectable class; `BEMECGElectro` and `pseudoECGElectro` are its concrete subclasses. `monoDomainElectro` gains `initialiseConductivityTensor()` (static utility, used by both monodomain and ecgModel subclasses) and `readVm()` (for postProcess mode).

**Tech Stack:** OpenFOAM C++ (wmake), existing runtime-selection pattern (`declareRunTimeSelectionTable`, `addToRunTimeSelectionTable`), `volTensorField`, `autoPtr`.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.H` | Add `initialiseConductivityTensor()` (protected static) and `readVm()` (protected) |
| Modify | `src/electroModels/monoDomainElectro/monoDomainElectro.C` | Implement both; refactor `initialiseConductivity()` to delegate |
| **Create** | `src/electroModels/ecgModel/ecgModel.H` | Abstract base + runtime-selector declaration |
| **Create** | `src/electroModels/ecgModel/ecgModel.C` | `New()` factory, `defineRunTimeSelectionTable` |
| **Create** | `src/electroModels/ecgElectro/ecgElectro.H` | Concrete electroModel: `postProcess_`, `ecgModelPtr_` |
| **Create** | `src/electroModels/ecgElectro/ecgElectro.C` | Constructor, `evolve()` mode switch |
| **Create** | `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.H` | Green's function ecgModel subclass |
| **Create** | `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.C` | `compute()`: Is_, Green's function integral, ECG.dat, VTK |
| **Create** | `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.H` | Gima-Rudy dipole ecgModel subclass |
| **Create** | `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.C` | `compute()`: dipole integral, pseudoECG.dat |
| Modify | `src/electroModels/Make/files` | Add all new classes; remove `greensFunctionECGElectro` |
| Delete | `src/electroModels/greensFunctionECGElectro/` | Replaced |
| Modify | `tutorials/ECG/constant/electroProperties` | Switch to `ecgElectro` + `ecgModel BEMECGElectro` |

`pdeECGElectro` is **not changed** in this plan — it stays as a direct `monoDomainElectro` subclass.

---

## Chunk 1: `monoDomainElectro` Infrastructure

**Goal:** Extract the conductivity tensor init pattern into a reusable static utility; add `readVm()` for postProcess mode. Existing behaviour unchanged — `NiedererEtAl2012` must still pass.

### Task 1: Add `initialiseConductivityTensor` and `readVm` to `monoDomainElectro`

**Files:**
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.H`
- Modify: `src/electroModels/monoDomainElectro/monoDomainElectro.C`

- [ ] **Step 1: Add declarations to `.H`**

  In the `protected` section, after the existing `conductivity()` accessor:

  ```cpp
  //- Initialise a conductivity tensor: try disk (READ_IF_PRESENT), fall back to dict.
  //  Used by both monoDomainElectro (for conductivity_) and ecgModel subclasses (for Gi_).
  static tmp<volTensorField> initialiseConductivityTensor
  (
      const word& fieldName,
      const dictionary& dict,
      const fvMesh& mesh
  );

  //- Re-read Vm from the current time directory.
  //  Used by ecgElectro in postProcess mode to skip the monodomain solve.
  void readVm();
  ```

- [ ] **Step 2: Implement `initialiseConductivityTensor` in `.C`**

  Add before the first constructor. Check the existing `initialiseConductivity()` for the exact dimensions used (look for `dimensionSet(-1,-3,3,0,0,2,0)` or similar) and use the same:

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
          fieldName, mesh.time().timeName(), mesh,
          IOobject::READ_IF_PRESENT, IOobject::NO_WRITE
      );

      // Start with a zero-valued field
      auto tresult = tmp<volTensorField>::New
      (
          IOobject
          (
              fieldName, mesh.time().timeName(), mesh,
              IOobject::NO_READ, IOobject::NO_WRITE
          ),
          mesh,
          dimensionedTensor(fieldName, io.headerOk()
              ? dimless  // will be overwritten below
              : dimensionedSymmTensor(fieldName, dict).dimensions(),
              tensor::zero)
      );

      if (io.headerOk())
      {
          tresult.ref() = volTensorField(io, mesh);
          Info<< "Read " << fieldName << " from disk" << nl;
      }
      else
      {
          const dimensionedSymmTensor s(fieldName, dict);
          tresult.ref() = dimensionedTensor(fieldName, s.dimensions(), s.value() & tensor::I);
          Info<< fieldName << " initialized from dict: " << s << nl;
      }

      return tresult;
  }
  ```

  Note: the above zero-field construction is illustrative — match the exact pattern from the existing `initialiseConductivity()` body. The goal is to extract that body verbatim into this static method and have `initialiseConductivity()` call it.

- [ ] **Step 3: Refactor `initialiseConductivity()` to delegate**

  Replace its body with:

  ```cpp
  tmp<volTensorField> monoDomainElectro::initialiseConductivity() const
  {
      return initialiseConductivityTensor("conductivity", cardiacProperties_, mesh());
  }
  ```

- [ ] **Step 4: Implement `readVm()`**

  ```cpp
  void Foam::monoDomainElectro::readVm()
  {
      Vm_.read();
  }
  ```

- [ ] **Step 5: Build and verify NiedererEtAl2012**

  Run `wmake` from the library root (triggers `wmakeLnInclude` to update symlinks):

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: no errors. Then:

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials/NiedererEtAl2012
  ./Allrun
  ```

  Expected: completes successfully; `postProcessing/` populated as before.

- [ ] **Step 6: Commit**

  ```bash
  git add src/electroModels/monoDomainElectro/
  git commit -m "refactor(monoDomainElectro): extract initialiseConductivityTensor + add readVm()"
  ```

---

## Chunk 2: `ecgModel` Abstract Base

**Goal:** New abstract class with runtime-selection — the `ionicModel` equivalent for ECG approaches.

### Task 2: Create `ecgModel`

**Files:**
- Create: `src/electroModels/ecgModel/ecgModel.H`
- Create: `src/electroModels/ecgModel/ecgModel.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `ecgModel.H`**

  Model this closely on `ionicModel.H`. Key elements:

  ```cpp
  #ifndef ecgModel_H
  #define ecgModel_H

  #include "monoDomainElectro.H"   // for initialiseConductivityTensor
  #include "volFields.H"
  #include "dimensionedScalar.H"
  #include "autoPtr.H"
  #include "runTimeSelectionTables.H"

  namespace Foam
  {

  class ecgModel
  {
  protected:

      const volScalarField& Vm_;    // Vm owned by ecgElectro (monoDomainElectro)
      const fvMesh&         mesh_;

      volTensorField        Gi_;    // intracellular conductivity (time-invariant)
      dimensionedScalar     sigmaT_; // torso conductivity (constant)

      // Electrode data common to both BEM and pseudo approaches
      wordList          electrodeNames_;
      List<vector>      electrodePositions_;

      void readElectrodes(const dictionary& dict);


  public:

      TypeName("ecgModel");

      declareRunTimeSelectionTable
      (
          autoPtr,
          ecgModel,
          dictionary,
          (
              const volScalarField& Vm,
              const dictionary& dict
          ),
          (Vm, dict)
      );

      static autoPtr<ecgModel> New
      (
          const volScalarField& Vm,
          const dictionary& dict
      );

      // Constructor
      ecgModel(const volScalarField& Vm, const dictionary& dict);

      virtual ~ecgModel() = default;

      // Access
      const volTensorField& Gi()     const { return Gi_; }
      const dimensionedScalar& sigmaT() const { return sigmaT_; }

      //- Compute ECG output for current timestep
      virtual void compute() = 0;

      virtual bool read(const dictionary& dict);
  };

  } // End namespace Foam

  #endif
  ```

  Note: `electrodeNames_` and `electrodePositions_` are placed in the base because both `BEMECGElectro` and `pseudoECGElectro` need electrodes. `readElectrodes()` reads the `electrodes` subdict verbatim from the existing `greensFunctionECGElectro` implementation.

- [ ] **Step 2: Write `ecgModel.C`**

  ```cpp
  #include "ecgModel.H"
  #include "addToRunTimeSelectionTable.H"

  namespace Foam
  {
      defineTypeNameAndDebug(ecgModel, 0);
      defineRunTimeSelectionTable(ecgModel, dictionary);
  }

  Foam::ecgModel::ecgModel
  (
      const volScalarField& Vm,
      const dictionary& dict
  )
  :
      Vm_(Vm),
      mesh_(Vm.mesh()),
      Gi_
      (
          monoDomainElectro::initialiseConductivityTensor("Gi", dict, Vm.mesh())
      ),
      sigmaT_("sigmaT", dict)
  {
      readElectrodes(dict);

      Info<< "ecgModel: Gi read, sigmaT = " << sigmaT_.value()
          << ", " << electrodeNames_.size() << " electrode(s)" << nl;
  }

  void Foam::ecgModel::readElectrodes(const dictionary& dict)
  {
      // Copy verbatim from greensFunctionECGElectro::readElectrodes()
      // Reads dict.subDict("electrodes"), populates electrodeNames_ and electrodePositions_
      if (dict.found("electrodes"))
      {
          const dictionary& elDict = dict.subDict("electrodes");
          electrodeNames_ = elDict.toc();
          electrodePositions_.setSize(electrodeNames_.size());
          forAll(electrodeNames_, i)
              electrodePositions_[i] = elDict.lookup(electrodeNames_[i]);
      }
  }

  Foam::autoPtr<Foam::ecgModel> Foam::ecgModel::New
  (
      const volScalarField& Vm,
      const dictionary& dict
  )
  {
      const word modelType(dict.lookup("ecgModel"));
      Info<< "Selecting ecgModel " << modelType << nl;

      auto* ctorPtr = dictionaryConstructorTable(modelType);
      if (!ctorPtr)
      {
          FatalErrorInFunction
              << "Unknown ecgModel type " << modelType << nl
              << "Valid types:" << dictionaryConstructorTablePtr_->sortedToc()
              << exit(FatalError);
      }

      return autoPtr<ecgModel>(ctorPtr(Vm, dict.subDict(modelType + "Coeffs")));
  }

  bool Foam::ecgModel::read(const dictionary& dict)
  {
      return true;
  }
  ```

- [ ] **Step 3: Add to `Make/files`**

  ```
  ecgModel/ecgModel.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly. No concrete subclasses yet.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/ecgModel/ src/electroModels/Make/files
  git commit -m "feat: add ecgModel abstract base with runtime-selection (ionicModel pattern)"
  ```

---

## Chunk 3: `ecgElectro` Concrete ElectroModel

**Goal:** The single user-facing `electroModel` for ECG. Owns the `ecgModel` autoPtr and the mode switch.

### Task 3: Create `ecgElectro`

**Files:**
- Create: `src/electroModels/ecgElectro/ecgElectro.H`
- Create: `src/electroModels/ecgElectro/ecgElectro.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `ecgElectro.H`**

  ```cpp
  #ifndef ecgElectro_H
  #define ecgElectro_H

  #include "monoDomainElectro.H"
  #include "ecgModel.H"

  namespace Foam
  {

  class ecgElectro
  :
      public monoDomainElectro
  {
      const Switch      postProcess_;
      autoPtr<ecgModel> ecgModelPtr_;

  public:

      TypeName("ecgElectro");

      ecgElectro
      (
          Time& runTime,
          const word& region = dynamicFvMesh::defaultRegion
      );

      virtual ~ecgElectro() = default;

      virtual bool evolve();
      virtual bool read();
  };

  } // End namespace Foam

  #endif
  ```

- [ ] **Step 2: Write `ecgElectro.C`**

  ```cpp
  #include "ecgElectro.H"
  #include "addToRunTimeSelectionTable.H"

  namespace Foam
  {
      defineTypeNameAndDebug(ecgElectro, 0);
      addToRunTimeSelectionTable(electroModel, ecgElectro, dictionary);
  }

  Foam::ecgElectro::ecgElectro
  (
      Time& runTime,
      const word& region
  )
  :
      monoDomainElectro(typeName, runTime, region),
      postProcess_
      (
          electroProperties().subDict(typeName + "Coeffs")
              .lookupOrDefault<Switch>("postProcess", false)
      ),
      ecgModelPtr_
      (
          ecgModel::New
          (
              Vm(),   // protected accessor from monoDomainElectro
              electroProperties().subDict(typeName + "Coeffs")
          )
      )
  {
      Info<< "ecgElectro: postProcess = " << postProcess_ << nl;
  }

  bool Foam::ecgElectro::evolve()
  {
      bool ok = true;

      if (!postProcess_)
          ok = monoDomainElectro::evolve();
      else
          readVm();

      ecgModelPtr_->compute();

      return ok;
  }

  bool Foam::ecgElectro::read()
  {
      return monoDomainElectro::read();
  }
  ```

  Note: `Vm()` — verify the exact name of the protected accessor that exposes `Vm_` in `monoDomainElectro`. If it doesn't exist as a non-const reference, add one or store a reference to the field before constructing the ecgModel.

- [ ] **Step 3: Add to `Make/files`**

  ```
  ecgElectro/ecgElectro.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: compiles cleanly. `ecgElectro` is registered in the runtime selector.

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/ecgElectro/ src/electroModels/Make/files
  git commit -m "feat: add ecgElectro concrete electroModel owning autoPtr<ecgModel>"
  ```

---

## Chunk 4: `BEMECGElectro` and `pseudoECGElectro` ecgModel Subclasses

**Goal:** Two concrete `ecgModel` subclasses extracted from `greensFunctionECGElectro`. Tutorial switched to `ecgElectro` + `BEMECGElectro` and verified end-to-end.

### Task 4: Create `BEMECGElectro`

**Files:**
- Create: `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.H`
- Create: `src/electroModels/ecgModels/BEMECGElectro/BEMECGElectro.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `BEMECGElectro.H`**

  ```cpp
  #ifndef BEMECGElectro_H
  #define BEMECGElectro_H

  #include "ecgModel.H"
  #include "OFstream.H"
  #include "triSurface.H"

  namespace Foam
  {

  class BEMECGElectro : public ecgModel
  {
      volScalarField        Is_;           // current source density

      autoPtr<triSurface>   torsoSurfacePtr_;
      pointField            torsoFaceCentres_;
      fileName              torsoVtkDir_;

      autoPtr<OFstream>     outputPtr_;    // ECG.dat

      void loadTorsoSurface(const dictionary& dict);
      void writeHeader();
      void writeTorsoVtk(const scalarList& phiGreens) const;

  public:

      TypeName("BEMECGElectro");

      BEMECGElectro(const volScalarField& Vm, const dictionary& dict);
      virtual ~BEMECGElectro() = default;

      virtual void compute();
  };

  } // End namespace Foam
  #endif
  ```

- [ ] **Step 2: Write `BEMECGElectro.C`**

  Constructor: copy from `greensFunctionECGElectro` constructor, removing `Gi_`/`sigmaT_` init (now in base) and pseudo-ECG stream setup. Call `loadTorsoSurface(dict)` if `dict.found("torsoSurface")`.

  `compute()`: copy the ECG block from `greensFunctionECGElectro::evolve()`:
  - Access Vm via `Vm_` (base member)
  - Access Gi via `Gi()` (base accessor)
  - Access sigmaT via `sigmaT()` (base accessor)
  - Compute `Is_ = -fvc::div(Gi() & fvc::grad(Vm_))`
  - Run **Green's function integral only** (drop pseudo-ECG integral entirely)
  - Write `ECG.dat` gated on `runTime().outputTime()` — access time via `mesh_.time()`
  - Write torso VTK if torso surface loaded

  `loadTorsoSurface`, `writeHeader`, `writeTorsoVtk`: copy verbatim from `greensFunctionECGElectro.C`, adapting member names.

  ```cpp
  namespace Foam
  {
      defineTypeNameAndDebug(BEMECGElectro, 0);
      addToRunTimeSelectionTable(ecgModel, BEMECGElectro, dictionary);
  }
  ```

- [ ] **Step 3: Add to `Make/files`**

  ```
  ecgModels/BEMECGElectro/BEMECGElectro.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/ecgModels/BEMECGElectro/ src/electroModels/Make/files
  git commit -m "feat: add BEMECGElectro ecgModel (Green's function ECG)"
  ```

---

### Task 5: Create `pseudoECGElectro`

**Files:**
- Create: `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.H`
- Create: `src/electroModels/ecgModels/pseudoECGElectro/pseudoECGElectro.C`
- Modify: `src/electroModels/Make/files`

- [ ] **Step 1: Write `pseudoECGElectro.H`**

  ```cpp
  #ifndef pseudoECGElectro_H
  #define pseudoECGElectro_H

  #include "ecgModel.H"
  #include "OFstream.H"

  namespace Foam
  {

  class pseudoECGElectro : public ecgModel
  {
      autoPtr<OFstream> outputPtr_;    // pseudoECG.dat

      void writeHeader();

  public:

      TypeName("pseudoECGElectro");

      pseudoECGElectro(const volScalarField& Vm, const dictionary& dict);
      virtual ~pseudoECGElectro() = default;

      virtual void compute();
  };

  } // End namespace Foam
  #endif
  ```

- [ ] **Step 2: Write `pseudoECGElectro.C`**

  Constructor: call `ecgModel(Vm, dict)` base (reads Gi_, sigmaT_, electrodes). Open `postProcessing/pseudoECG.dat`. Call `writeHeader()`.

  `compute()`:

  ```cpp
  void Foam::pseudoECGElectro::compute()
  {
      const auto tgradVm = fvc::grad(Vm_);
      const auto& gradVm = tgradVm();

      const scalarField& Vols = mesh_.V();
      const vectorField& C    = mesh_.C();
      const label nPts        = electrodePositions_.size();

      scalarList localSums(nPts, 0.0);

      forAll(C, cI)
      {
          forAll(electrodePositions_, pI)
          {
              const vector r_vec = C[cI] - electrodePositions_[pI];
              const scalar r = Foam::mag(r_vec) + SMALL;
              const vector dipole = (Gi()[cI] & gradVm[cI]) * Vols[cI];
              localSums[pI] -= (dipole & r_vec) / (r*r*r);
          }
      }

      Foam::reduce(localSums, sumOp<scalarList>());

      if (mesh_.time().outputTime() && Pstream::master())
      {
          *outputPtr_ << mesh_.time().timeOutputValue();
          for (const scalar v : localSums)
              *outputPtr_ << tab << v;
          *outputPtr_ << nl;
      }
  }
  ```

  Check existing `greensFunctionECGElectro.C` ~lines 360–390 for the exact dipole formula and match it precisely (including any `sigmaT` factor or lack thereof).

  ```cpp
  namespace Foam
  {
      defineTypeNameAndDebug(pseudoECGElectro, 0);
      addToRunTimeSelectionTable(ecgModel, pseudoECGElectro, dictionary);
  }
  ```

- [ ] **Step 3: Add to `Make/files`**

  ```
  ecgModels/pseudoECGElectro/pseudoECGElectro.C
  ```

- [ ] **Step 4: Build**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/electroModels/ecgModels/pseudoECGElectro/ src/electroModels/Make/files
  git commit -m "feat: add pseudoECGElectro ecgModel (Gima-Rudy dipole)"
  ```

---

### Task 6: Update tutorial and verify end-to-end

**Files:**
- Modify: `tutorials/ECG/constant/electroProperties`

- [ ] **Step 1: Update `electroProperties`**

  Full replacement:

  ```
  electroModel  ecgElectro;

  ecgElectroCoeffs
  {
      postProcess   false;

      // Monodomain settings (verbatim from greensFunctionECGElectroCoeffs):
      chi               [0 -1 0 0 0 0 0]    140000;
      Cm                [-1 -4 4 0 0 2 0]   0.01;
      ionicModel        BuenoOrovio;
      tissue            epicardialCells;
      solutionAlgorithm explicit;
      solver            RKF45;
      initialODEStep    1e-6;
      maxSteps          1000000000;

      monodomainStimulus
      {
          stimulusLocationMin  (0 0 5.5e-3);
          stimulusLocationMax  (1.5e-3 1.5e-3 7e-3);
          stimulusDuration     [0 0 1 0 0 0 0] 2e-3;
          stimulusIntensity    [0 -3 0 0 0 1 0] 50000;
          stimulusStartTime    0.0;
      }

      // ECG type selection:
      ecgModel  BEMECGElectro;

      BEMECGElectroCoeffs
      {
          Gi      [-1 -3 3 0 0 2 0] (0.17 0 0  0.019 0  0.019);
          sigmaT  [-1 -3 3 0 0 2 0] 0.24725;

          electrodes
          {
              V1 (-26.1184 -283.543 -70.0558);
              V2 ( 73.582  -278.089 -70.1427);
              V3 ( 85.1153 -276.414 -93.8405);
              V4 ( 99.7976 -269.695 -106.596);
              V5 (122.885  -253.215 -119.155);
              V6 (148.613  -227.84  -119.317);
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
  - `postProcessing/ECG.dat` populated with time series for V1–V6

- [ ] **Step 3: Commit**

  ```bash
  git add tutorials/ECG/constant/electroProperties
  git commit -m "feat(tutorials/ECG): switch to ecgElectro + BEMECGElectro"
  ```

---

## Chunk 5: Cleanup

**Goal:** Remove `greensFunctionECGElectro`, run full regression.

### Task 7: Remove `greensFunctionECGElectro` and final regression

- [ ] **Step 1: Remove from `Make/files`**

  Delete the line:
  ```
  greensFunctionECGElectro/greensFunctionECGElectro.C
  ```

- [ ] **Step 2: Delete directory**

  ```bash
  rm -rf /Users/simaocastro/cardiacFoamv2/src/electroModels/greensFunctionECGElectro
  ```

- [ ] **Step 3: Build clean**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/src/electroModels
  wmake
  ```

  Expected: no references to `greensFunctionECGElectro` anywhere.

- [ ] **Step 4: Full regression**

  ```bash
  cd /Users/simaocastro/cardiacFoamv2/tutorials
  ./Alltest-regression
  ```

  Expected: all tutorials pass.

- [ ] **Step 5: Commit**

  ```bash
  git add -A src/electroModels/
  git commit -m "chore: remove greensFunctionECGElectro (replaced by ecgElectro + BEMECGElectro)"
  ```

---

## Deferred / Known Limitations

1. **`Gi` tensor → scalar simplification**: `Gi_` stays as a full anisotropic `volTensorField`. The relationship `Gi = sigmaI * conductivity()` is noted but deferred.

2. **`postProcess` time-looping**: `readVm()` is wired up in `ecgElectro::evolve()`, but the solver binary's time-directory looping for standalone postProcess invocation is not implemented in this plan.

3. **`pdeECGElectro` not migrated**: stays as a direct `monoDomainElectro` subclass with its own `Gi_`/`Ge_`/`sigmaT_`. Migration to use `ecgModel` is deferred due to the separate torso mesh complexity.
