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
    manufacturedGraphVerifier

Description
    Manufactured monodomain verifier for the 1D graph.
    Calculates L1, L2, and Linf errors against exact analytical
    solutions across the 1D graph nodes.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedGraphVerifier_H
#define manufacturedGraphVerifier_H

#include "graphVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedGraphVerifier : public graphVerificationModel
{
    mutable bool errorsReported_;

public:

    OverrideTypeName("manufacturedGraphVerifier");

    explicit manufacturedGraphVerifier(const dictionary& dict);

    virtual void preProcess
    (
        const Time& runTime,
        ionicModel* model,
        scalarField& Vm1D,
        const pointField& nodeLocations
    ) override;

    virtual void postProcess
    (
        const Time& runTime,
        const ionicModel* model,
        const scalarField& Vm1D,
        const pointField& nodeLocations
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
