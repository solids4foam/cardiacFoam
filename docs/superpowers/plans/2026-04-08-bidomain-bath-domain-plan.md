# Bidomain-Bath Domain Audit And Implementation Plan

## Goal

Add a true bidomain-with-bath workflow that follows the existing domain/solver/coupler architecture instead of stretching the current pseudo-ECG path beyond its contract.

This note also records one preparatory code change already made:

- `electroStateProvider` and `myocardiumSolver` now expose

  `intracellularConductivityPtr()` and `extracellularConductivityPtr()`.

- `BidomainSolver` maps those to `Gi` and `Ge`.

- Existing `conductivityPtr()` behavior is unchanged, so the current

  pseudo-ECG path still sees the same tensor it used before.

## Audit Summary

### What exists today

- Primary tissue domain:

  `src/electroModels/electroDomains/myocardiumDomain/`

- Downstream scalar ECG domain:

  `src/electroModels/electroDomains/ecgDomain/`

- ECG kernels:

  `src/electroModels/ecgModels/`

- Generic one-way state sharing:

  `src/electroModels/core/electroStateProvider.H`

- Generic domain-to-domain coupling hooks:

  `src/electroModels/electroCouplers/`

### Why the current ECG path is not enough

The current `ECGDomain` is designed for downstream post-processing on the
upstream mesh:

- it binds its mesh to `stateProvider.mesh()`

- it owns scalar output buffers, not a bath potential field

- its solver API is `compute(domain, scalarField&)`

- it reads upstream state directly instead of receiving prepared interface data

That is appropriate for `pseudoECG`, but not for a bath PDE solve on a
separate bath mesh.

### What the FDA document requires

The bath model is a genuine PDE in the bath region `Omega_b`:

- in tissue `Omega`, the bidomain equations hold

- in bath `Omega_b`, `div(sigma_b grad(phi_e)) = 0`

- at the heart-bath interface:

  - continuity of `phi_e`
  - continuity of extracellular current

- at the outer bath boundary:

  - either Neumann surface stimulus
  - or a Dirichlet ground electrode

This means the implementation needs:

- a bath-owned field for `phi_e` on the bath mesh

- a heart-to-bath interface mapping step

- access to extracellular conductivity `Ge` on the heart side

## Recommended Architecture

### 1. New downstream domain

Add:

- `src/electroModels/electroDomains/bathDomain/bathDomain.H`

- `src/electroModels/electroDomains/bathDomain/bathDomain.C`

- optional `src/electroModels/electroDomains/bathDomain/README.md`

Responsibilities:

- own the bath mesh-backed state

- own the bath solution field, e.g. `phiEbath_`

- own bath material properties, e.g. `sigmaBath_`

- implement `electroDomainInterface`

- optionally implement `electroStateProvider` if later domains need to read

  bath potentials

- delegate the numerical kernel to a runtime-selected solver

Suggested class:

- `BathDomain`

Suggested key dictionary block:

- `bathDomains { <name> { ... } }`

### 2. Keep the bath PDE kernel under `ecgModels`

For the current phase, do not introduce a new top-level `bathModels/`
directory. The bath-region PDE is still part of the body-potential / ECG side
of the codebase, so keep the kernel under `ecgModels`.

Add:

- `src/electroModels/ecgModels/bidomainBathECGSolver/bidomainBathECGSolver.H`

- `src/electroModels/ecgModels/bidomainBathECGSolver/bidomainBathECGSolver.C`

Responsibilities:

- `bidomainBathECGSolver` solves the bath-region elliptic problem for `phi_e`

- the owning `BathDomain` still remains separate and owns the bath-side fields

  and lifecycle

Suggested runtime type:

- `bidomainBathECG`

Important design note:

- keep `BathDomain` as a separate domain now

- keep `bidomainBathECGSolver` under `ecgModels/`

- postpone any broader refactor of the current `ECGDomain` until later

That gives a minimal first step:

- new domain for the bath region

- ECG-side numerical kernel for body-potential computation

- no forced redesign of the current pseudo-ECG path

### 3. New bath interface endpoint

Add to:

- `src/electroModels/electroCouplers/electroDomainCouplingEndpoints.H`

New endpoint idea:

```cpp
class bathCouplingEndpoint
{
public:
    virtual ~bathCouplingEndpoint() = default;

    virtual const fvMesh& mesh() const = 0;

    virtual volScalarField& bathPotentialField() = 0;

    virtual void setInterfaceSource(const scalarField& values) = 0;
};

```

Exact shape depends on the discretization choice:

- mapped Dirichlet data on the bath-side interface

- mapped Neumann flux data on the bath-side interface

- or an assembled interface source object

Do not encode those details in `electroStateProvider`; that interface is
read-only state exposure, not a coupling contract.

### 4. New downstream coupler

Add:

- `src/electroModels/electroCouplers/heartBathInterfaceCoupler.H`

- `src/electroModels/electroCouplers/heartBathInterfaceCoupler.C`

Responsibilities:

- read the updated heart-side fields after the myocardium solve

- pull `phiE` and extracellular conductivity from the myocardium provider

- map interface data from the heart boundary to the bath boundary

- prepare the coupling consumed by `BathDomain`

This coupler should use `preparePostPrimaryCoupling()`, because the bath is a
downstream domain advanced after the primary myocardium domain.

## Why `electroStateProvider` Is Useful But Not Sufficient

`electroStateProvider` should still be used for:

- `VmPtr()`

- `phiEPtr()`

- `intracellularConductivityPtr()`

- `extracellularConductivityPtr()`

That gives downstream domains and couplers read-only access to heart-side
fields.

It is not sufficient by itself because it does not:

- own bath-side state

- represent a separate bath mesh

- carry interface mapping logic

- prepare boundary data between disjoint domains

That work belongs in a dedicated downstream coupler.

## Concrete File Changes To Implement

### Core interfaces

Already changed:

- `src/electroModels/core/electroStateProvider.H`

- `src/electroModels/electroDomains/myocardiumDomain/myocardiumSolver.H`

- `src/electroModels/electroDomains/myocardiumDomain/myocardiumDomain.H`

- `src/electroModels/myocardiumModels/bidomainSolver/bidomainSolver.H`

Future likely changes:

- `src/electroModels/electroCouplers/electroDomainCouplingEndpoints.H`

### New domain

- `src/electroModels/electroDomains/bathDomain/bathDomain.H`

- `src/electroModels/electroDomains/bathDomain/bathDomain.C`

### New ECG-side bath solver

- `src/electroModels/ecgModels/bidomainBathECGSolver/bidomainBathECGSolver.H`

- `src/electroModels/ecgModels/bidomainBathECGSolver/bidomainBathECGSolver.C`

### New coupler

- `src/electroModels/electroCouplers/heartBathInterfaceCoupler.H`

- `src/electroModels/electroCouplers/heartBathInterfaceCoupler.C`

### Builder/orchestration integration

- `src/electroModels/core/electrophysicsSystemBuilder.H`

- `src/electroModels/core/electrophysicsSystemBuilder.C`

- `src/electroModels/core/electrophysicsSystem.H`

- `src/electroModels/core/electrophysicsSystem.C`

The existing downstream-domain list can likely be reused, but the builder needs
new discovery/configuration helpers for bath domains and bath couplers.

### Build registration

- `src/electroModels/Make/files`

- possibly `src/electroModels/Make/options`

Add new compilation units for:

- `bathDomain`

- `bidomainBathECGSolver`

- `heartBathInterfaceCoupler`

## Suggested Dictionary Shape

```foam
bathDomains
{
    torsoBath
    {
        bathSolver  bidomainBathECG;

        sigmaBath   [ ... ] 0.5;

        outerBoundaryCondition  grounded;
        groundPatches (leftElectrode);

        electrodeStimulus
        {
            mode surfaceCurrent;
            patches (leftElectrode rightElectrode);
            values  (-0.01 0.01);
        }
    }
}

bathCouplings
{
    heartToBath
    {
        electroDomainCoupler  heartBathInterfaceCoupler;
        primaryDomain myocardium;
        bathDomain torsoBath;
        couplingMode continuityPhiEAndFlux;
    }
}

```

This follows the existing builder pattern:

- domain collections are declared by named subdictionaries

- couplers are declared separately by named subdictionaries

## Minimal Implementation Sequence

1. Add `BathDomain`.
2. Add `bidomainBathECGSolver` under `ecgModels/` so it solves the bath
   elliptic PDE with simple
   outer BC support.
3. Add `bathCouplingEndpoint` and `heartBathInterfaceCoupler`.
4. Extend the builder to read `bathDomains` and `bathCouplings`.
5. Add a manufactured tutorial case that mirrors the FDA 1D setup first.

Start with the 1D manufactured-with-bath problem only. The FDA document notes
the 2D and 3D bath cases are still effectively 1D in the solution structure,
so a 1D-first implementation is the right lowest-risk path.

## Design Constraints To Preserve

- Keep domain ownership separate from solver kernels.

- Keep interface mapping in couplers, not in domains.

- Keep read-only field sharing in `electroStateProvider`.

- For now, do not overload `ECGDomain` with a separate-mesh bath PDE

  responsibility.

- It is acceptable for the bath PDE kernel itself to live under `ecgModels/`.

- Preserve the existing `pseudoECG` path unchanged.
