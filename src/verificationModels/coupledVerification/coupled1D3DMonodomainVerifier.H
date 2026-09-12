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

Class
    coupled1D3DMonodomainVerifier

Description
    Diagnostic verifier for the staged 1D-3D manufactured monodomain workflow.

    This verifier operates at the coupler level, meaning it receives both
    the primary domain (3D myocardium) and secondary domain (1D graph) to
    compare their values against the analytical coupled solutions.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef coupled1D3DMonodomainVerifier_H
#define coupled1D3DMonodomainVerifier_H

#include "couplingVerificationModel.H"
#include "overrideTypeName.H"
#include "scalarField.H"

namespace Foam
{

class coupled1D3DMonodomainVerifier : public couplingVerificationModel
{
    mutable bool diagnosticsWritten_;
    mutable scalarField exactPrimarySource_;
    mutable scalarField exactSecondaryAppliedCurrent_;

public:

    OverrideTypeName("coupled1D3DMonodomainVerifier");

    explicit coupled1D3DMonodomainVerifier(const dictionary& dict);

    virtual void preProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) override;

    virtual void updateManufacturedSource
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        scalar primaryTime,
        scalar secondaryTime,
        bool implicitCoupling,
        bool bidirectionalCoupling,
        const word& phaseName
    ) override;

    virtual void postProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
