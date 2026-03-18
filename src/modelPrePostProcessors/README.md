# modelPrePostProcessors library architecture

This directory provides reusable model pre/post processors that can be shared
across electro, active-tension, and ECG model stacks.

## Directory structure

```text
src/modelPrePostProcessors/
├── electroModelsPrePostProcessor/      # Generic electro-model pre/post interface
├── tmanufacturedFDAPrePostProcessor/   # Manufactured-solution verification processor
├── Make/
└── README.md
```

## Design intent

- Keep reusable pre/post lifecycle abstractions outside any one model family.
- Allow electro-model processors today, while leaving room for future
  active-tension or ECG-specific processors in the same library.
- Keep analytical manufactured-solution helpers under `ionicModels`, while
  solver-facing verification/reporting code lives here.
