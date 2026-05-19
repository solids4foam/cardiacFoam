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
    electroDomainCoupler

Description
    Shared context for coupling models that exchange data between domains.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef electroDomainCoupler_H
#define electroDomainCoupler_H

#include "autoPtr.H"
#include "electroDomainCouplingEndpoints.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class electroDomainInterface;

class electroDomainCoupler
{
protected:

    const fvMesh& mesh_;
    tissueCouplingEndpoint& primaryDomain_;
    electroDomainInterface& secondaryDomain_;

public:

    TypeName("electroDomainCoupler");

    declareRunTimeSelectionTable
    (
        autoPtr,
        electroDomainCoupler,
        dictionary,
        (
            tissueCouplingEndpoint& primaryDomain,
            electroDomainInterface&       secondaryDomain,
            const dictionary&         dict
        ),
        (primaryDomain, secondaryDomain, dict)
    );

    electroDomainCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface&       secondaryDomain
    )
    :
        mesh_(primaryDomain.mesh()),
        primaryDomain_(primaryDomain),
        secondaryDomain_(secondaryDomain)
    {}

    virtual ~electroDomainCoupler() = default;

    static autoPtr<electroDomainCoupler> New
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface&       secondaryDomain,
        const dictionary&         dict
    );

    //- Prepare coupling terms consumed by the secondary advance.
    //  This hook runs before the auxiliary/secondary domain is advanced.
    virtual void prepareSecondaryCoupling
    (
        scalar t0,
        scalar dt
    )
    {
        (void)t0;
        (void)dt;
    }

    //- Prepare coupling terms consumed by the primary myocardium advance.
    //  This hook runs after the network-domain advance and before the
    //  myocardium reaction/diffusion solve. The default keeps backward
    //  compatibility with older single-phase coupling models.
    virtual void preparePrimaryCoupling
    (
        scalar t0,
        scalar dt
    )
    {
        prepareCoupling(t0, dt);
    }

    //- Legacy single-phase coupling hook.
    virtual void prepareCoupling
    (
        scalar t0,
        scalar dt
    )
    {
        (void)t0;
        (void)dt;
    }

    //- Prepare coupling terms consumed by a post-primary domain.
    //  This hook runs after the myocardium solve and before a one-way
    //  post-primary auxiliary domain is advanced.
    virtual void preparePostPrimaryCoupling
    (
        scalar t0,
        scalar dt
    )
    {
        (void)t0;
        (void)dt;
    }

    //- Optional write hook.
    virtual void write() {}

    //- Optional teardown hook.
    virtual void end() {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
