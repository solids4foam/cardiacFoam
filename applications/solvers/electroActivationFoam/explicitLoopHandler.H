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

#ifndef explicitLoopHandler_H
#define explicitLoopHandler_H

#include "ionicModel.H"


using namespace Foam;


class explicitLoopHandler
{
private:
    ionicModel& ionicModel_;

    // Cached solver parameters
    scalar dx_;
    int dim_;
    scalar Def_max_;

    // Internal helper: compute dt stability limit
    scalar computeStableDt(const scalar CFL) const;

public:
    explicitLoopHandler(const fvMesh& mesh, ionicModel& model);

    // Precompute dx, dim, Def_max and update dt
    void initializeExplicit
    (
        scalar& dt,
        int& nsteps,
        const scalar chiVal,
        const scalar CmVal,
        const volTensorField& conductivity,
        const scalar CFL,
        const scalar dx,
        const int dim

    );

    // Perform one explicit time step
    void explicitLoop
    (
        const scalar t0,
        const scalar dt,
        volScalarField& Vm,
        volScalarField& Iion,
        Field<Field<scalar>>& states,
        volScalarField& externalStimulusCurrent,
        const labelList& stimulusCellIDs,
        const scalar stimulusIntensity,
        const scalar stimulusDuration,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const volTensorField& conductivity
    );
};

#endif
