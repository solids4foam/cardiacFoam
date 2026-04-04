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

#ifndef myocardialDomainModel_H
#define myocardialDomainModel_H

#include "electroModel/electrophysicsDomain.H"
#include "electroModel/electroStateProvider.H"
#include "electromechanicalSignalProvider.H"
#include "electroDomainCouplingEndpoints.H"
#include "pimpleControl.H"

namespace Foam
{

class myocardialDomainModel
:
    public electrophysicsDomain,
    public tissueCouplingEndpoint,
    public electroStateProvider
{
protected:

    const fvMesh& mesh_;

public:

    explicit myocardialDomainModel(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    virtual ~myocardialDomainModel() = default;

    virtual const fvMesh& mesh() const override
    {
        return mesh_;
    }

    virtual const Time& time() const override
    {
        return mesh_.time();
    }

    //- Field receiving system-level source terms for this domain.
    virtual volScalarField& sourceField() override = 0;

    virtual const volScalarField& Vm() const override = 0;

    //- Access to optional ionic/coupling signals produced by this domain.
    virtual const ElectromechanicalSignalProvider* provider() const = 0;

    using electrophysicsDomain::advance;

    //- Single integration step for the myocardium.
    //  When a pimpleControl pointer is provided, the domain should use
    //  implicit diffusion, otherwise explicit.
    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) = 0;

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const = 0;

    virtual bool shouldPostProcess() const = 0;

    virtual void exportStates() = 0;

    virtual void exportPostProcessFields() = 0;

    virtual void postProcess() = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
