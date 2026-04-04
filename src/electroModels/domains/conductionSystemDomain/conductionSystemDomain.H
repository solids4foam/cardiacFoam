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
    Foam::ConductionSystemDomain

Description
    Abstract base class for conduction-network electro domains that expose
    terminal coupling endpoints to the primary myocardium.

    Concrete subclasses expose terminal endpoints and any domain-specific
    state required by the selected conduction-system solver. Coupling exchange
    with the primary tissue domain is handled by a separate domain-coupling
    model.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ConductionSystemDomain.C
\*---------------------------------------------------------------------------*/

#ifndef ConductionSystemDomain_H
#define ConductionSystemDomain_H

#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "electroModel/electrophysicsDomain.H"
#include "conductionSystemDomainContext.H"
#include "electroDomainCouplingEndpoints.H"

namespace Foam
{

class ConductionSystemDomain
:
    public electrophysicsDomain,
    public networkCouplingEndpoint
{
protected:

    ConductionSystemDomainContext context_;

public:

    TypeName("ConductionSystemDomain");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ConductionSystemDomain,
        dictionary,
        (
            const ConductionSystemDomainContext& context,
            const dictionary&     dict,
            const scalar          initialDeltaT
        ),
        (context, dict, initialDeltaT)
    );

    ConductionSystemDomain
    (
        const ConductionSystemDomainContext& context,
        const dictionary&     dict
    );

    static autoPtr<ConductionSystemDomain> New
    (
        const ConductionSystemDomainContext& context,
        const dictionary&                    dict,
        scalar                               initialDeltaT
    );

    virtual ~ConductionSystemDomain() = default;

    virtual void advance
    (
        scalar          t0,
        scalar          dt
    ) override = 0;

    virtual const Time& time() const override
    {
        return context_.time();
    }

    const ConductionSystemDomainContext& context() const
    {
        return context_;
    }

    virtual const labelList& terminalNodes() const override = 0;

    virtual const pointField& terminalLocations() const override = 0;

    virtual void terminalVm(scalarField& values) const override = 0;

    //- Receive terminal exchange terms prepared by a coupling model.
    //  Positive terminalCurrent removes charge from the network and injects
    //  it into the primary myocardium.
    virtual void setTerminalCoupling
    (
        const scalarField& terminalCurrent,
        const scalarField& terminalSource
    ) override = 0;

    virtual void write() override {}

    virtual void end() override {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
