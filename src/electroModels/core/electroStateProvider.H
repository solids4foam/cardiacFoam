/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef electroStateProvider_H
#define electroStateProvider_H

#include "volFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class electroStateProvider
                       ~~~~~~~~~~~~~~~~~~~~~~~

Description
    Abstract read-only state interface that enables one-way coupling between
    electrophysiology domains while maintaining clean separation of concerns.

    This interface is the cornerstone of the domain separation architecture:
    - **Post-primary domains** (e.g., ECG) read electrophysiology state
      without providing feedback to the primary domain (myocardium)
    - **Dependency inversion**: Domains depend on this abstract interface,
      not on concrete solver implementations (monodomain, bidomain, etc.)
    - **Optional state channels**: Via pointer-based opt-in pattern, domains
      gracefully handle missing fields (e.g., phiE unavailable in monodomain)

Design Pattern: Pointer-Based Opt-In
    Methods return const X* (pointer) instead of const X& (reference).
    Returning nullptr signals "this state channel is not available."
    Consumers check for nullptr before dereferencing.

    Example usage:
        const auto* VmPtr = stateProvider.VmPtr();
        if (VmPtr) {
            // Safe to use transmembrane voltage
            compute_ecg(*VmPtr);
        }

Benefits
    - Breaks cyclic dependencies between myocardium and ECG domains
    - Enables independent evolution schedules (ECG reads, doesn't block)
    - Improves testability (mock providers for unit tests)
    - Supports future solvers automatically (just inherit this interface)

Implementers
    Producers:   electroModel, MyocardiumDomain
    Consumers:   ECGDomain, ecgVerificationModel

See also
    ARCHITECTURE.md — Comprehensive domain separation design rationale

\*---------------------------------------------------------------------------*/

class electroStateProvider
{
public:

    virtual ~electroStateProvider() = default;

    //- The mesh on which this state provider is defined.
    virtual const fvMesh& mesh() const = 0;

    //- Optional state channels exposed by the provider.
    virtual const volScalarField* VmPtr() const
    {
        return nullptr;
    }

    virtual const volScalarField* phiEPtr() const
    {
        return nullptr;
    }

    virtual const volTensorField* conductivityPtr() const
    {
        return nullptr;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
