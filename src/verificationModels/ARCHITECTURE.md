# verificationModels

`src/verificationModels` holds six verifier families, each a concrete
implementation of an abstract, runtime-selected base class, plus the
manufactured/reference helpers they use.

## Directory layout

```text
src/verificationModels/
├── monodomainVerification/   # Monodomain manufactured/reference verifiers
├── bidomainVerification/     # Bidomain manufactured/reference verifiers
├── bathBidomainVerification/ # Bath-bidomain manufactured/reference verifiers
├── eikonalVerification/      # Eikonal manufactured/reference verifiers
├── ecgVerification/          # ECG verifier family (concrete verifiers only)
├── electromechanicsVerification/ # Electromechanics verifiers
├── coupledVerification/      # Coupled domain verifiers
├── Make/
└── README.md
```

## Verifier families

Every abstract base is a runtime-selection root, constructed from dictionary
input. Five of the six live in `electroModels/core/verificationModels/`; only
`electromechanicalVerificationModel` is defined locally, in this library's own
`electromechanicsVerification/` folder.

| Abstract base | Defined in | Concrete verifiers |
|---|---|---|
| `electroVerificationModel` | `electroModels/core/verificationModels/` | `manufacturedFDAMonodomainVerifier`, `manufacturedAnisotropicMonodomainVerifier`, `manufacturedFDABidomainVerifier`, `manufacturedFDABathBidomainVerifier` |
| `eikonalVerificationModel` | `electroModels/core/verificationModels/` | `manufacturedEikonalVerifier` |
| `ecgVerificationModel` | `electroModels/core/verificationModels/` | `manufacturedPseudoECGVerifier`, `manufacturedEikonalECGVerifier` |
| `couplingVerificationModel` | `electroModels/core/verificationModels/` | `coupled1D3DMonodomainVerifier` |
| `graphVerificationModel` | `electroModels/core/verificationModels/` | `manufacturedGraphVerifier` |
| `electromechanicalVerificationModel` | `verificationModels/electromechanicsVerification/` (local) | `manufacturedElectromechanicsVerifier` |

`electroVerificationModel` is the myocardium-side base: it provides
verification hooks that monodomain/bidomain spatial workflows call before or
after the main solve. `ecgVerificationModel` is built around an upstream
`electroStateProvider`. `couplingVerificationModel` covers bidirectional or
unidirectional coupled workflows. `electromechanicalVerificationModel` adds
`configured(dict)`, reporting whether a verification block is present.

## Manufactured/reference helpers

Shared analytical helpers — exact/reference fields and formulas used by the
concrete verifiers — sit beside each family's verifiers:

- `monodomainVerification/manufacturedFDAReference.H`
- `monodomainVerification/manufacturedAnisotropicMonodomainReference.H`
- `bidomainVerification/manufacturedFDABidomainReference.H`
- `bathBidomainVerification/manufacturedFDABathBidomainReference.H`
- `eikonalVerification/manufacturedEikonalReference.H`
- `electromechanicsVerification/manufacturedElectromechanicsReference.H`

## Notes on naming

The folder/runtime ionic-model name for a manufactured model
(`monodomainFDAManufactured`, etc., in `src/ionicModels`) and the
verification-model name for its verifier
(`manufacturedFDAMonodomainVerifier`, etc., here) are different layers —
don't conflate them.

## What this layer validates

Monodomain, bidomain, eikonal, and coupled electromechanics spatial
workflows; ECG post-processing (pseudo-ECG and eikonal-ECG); coupled 1D/3D
workflows; and 1D graph (conduction-system) solver workflows — through
manufactured/reference formulas, error reporting, and runtime-selected
verifier objects integrated into the main electro stack.

This folder does not own the myocardium or ECG solvers. It provides the
verification-side models called from those workflows.
