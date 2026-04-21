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

#ifndef reactionDiffusionPvjCoupler_H
#define reactionDiffusionPvjCoupler_H

#include "pvjCoupler.H"
#include "overrideTypeName.H"
#include "Switch.H"

namespace Foam
{

class ReactionDiffusionPvjCoupler
:
    public PVJCoupler
{
    scalar R_pvj_;
    Switch debugCoupling_;
    mutable scalarField tissueVmBuffer_;
    mutable scalarField networkVmBuffer_;

    void couplingCurrentAtPvjs
    (
        const scalarField& networkVm,
        const scalarField& tissueVm,
        scalarField& current
    ) const;

    void evaluateCoupling(const char* phaseName) const;

    void reportCouplingDiagnostics(const char* phaseName) const;

public:

    OverrideTypeName("reactionDiffusionPvjCoupler");

    ReactionDiffusionPvjCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        const dictionary& dict
    );

    virtual void prepareSecondaryCoupling(scalar t0, scalar dt) override;

    virtual void preparePrimaryCoupling(scalar t0, scalar dt) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
