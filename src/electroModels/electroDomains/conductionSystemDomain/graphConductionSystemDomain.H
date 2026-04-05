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

#ifndef GraphConductionSystemDomain_H
#define GraphConductionSystemDomain_H

#include "conductionSystemDomain.H"
#include "graphConductionSystemSolver.H"
#include "conductionGraph.H"
#include "overrideTypeName.H"

namespace Foam
{

class ionicModel;

class GraphConductionSystemDomain
:
    public ConductionSystemDomain
{
    autoPtr<GraphConductionSystemSolver> solverPtr_;

public:

    GraphConductionSystemDomain
    (
        const ConductionSystemDomainContext& context,
        const dictionary& dict
    );

    virtual ~GraphConductionSystemDomain() = default;

    virtual scalarField& membranePotential() = 0;

    virtual scalarField& ionicCurrent() = 0;

    virtual ionicModel& ionicModelRef() = 0;

    virtual const labelList& edgeStartNodes() const = 0;

    virtual const labelList& edgeEndNodes() const = 0;

    virtual const scalarField& edgeLengths() const = 0;

    virtual const scalarField& edgeConductances() const = 0;

    virtual scalar chi() const = 0;

    virtual scalar Cm() const = 0;

    virtual const conductionGraph& graph() const = 0;

    virtual scalarField& activationTime() = 0;

    virtual void assembleAppliedCurrent
    (
        scalar t0,
        scalarField& appliedCurrent
    ) const = 0;

    virtual void reportAdvanceDiagnostics
    (
        scalar t0,
        scalar dt
    ) const;

    virtual void advance
    (
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
