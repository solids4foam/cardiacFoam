# verificationModels

Built-in verification: verifiers that compare a run with a manufactured solution and report the error. They are runtime-selected like every other model, so verification runs in the normal solver.

## What's available

| Folder | Verifiers |
|---|---|
| `monodomainVerification/` | `manufacturedFDAMonodomainVerifier`, `manufacturedGraphVerifier`, `manufacturedAnisotropicMonodomainVerifier` |
| `bidomainVerification/` | `manufacturedFDABidomainVerifier` |
| `bathBidomainVerification/` | `manufacturedFDABathBidomainVerifier` |
| `eikonalVerification/` | `manufacturedEikonalVerifier` |
| `ecgVerification/` | `manufacturedPseudoECGVerifier`, `manufacturedEikonalECGVerifier`, `manufacturedLeadFieldECGVerifier` |
| `coupledVerification/` | `coupled1D3DMonodomainVerifier` |
| `electromechanicsVerification/` | `manufacturedElectromechanicsVerifier` |

Cases that run them are in [tutorials/manufacturedSolutions](../../tutorials/manufacturedSolutions/README.md). The monodomain and eikonal `*Reference.H` headers are generated from the sympy generators in [tutorials/manufacturedSolutions/generators](../../tutorials/manufacturedSolutions/generators/README.md); edit the generator, not the header ([overview](../../docs/manufactured-solutions.md)).

**Deep dive:** [ARCHITECTURE.md](ARCHITECTURE.md) explains how this library is built inside.

## What this does not own

- The abstract verifier bases. They live in `electroModels/core/verificationModels`, which is why this library builds after `electroModels`.
- The manufactured cell models: [ionicModels](../ionicModels/README.md).
