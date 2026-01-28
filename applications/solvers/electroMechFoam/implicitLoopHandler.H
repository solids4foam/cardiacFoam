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

#ifndef implicitLoopHandler_H
#define implicitLoopHandler_H

#include "ionicModel.H"
#include "pimpleControl.H"

using namespace Foam;


class implicitLoopHandler
{
private:
    ionicModel& ionicModel_;

public:
    explicit implicitLoopHandler(const fvMesh& mesh, ionicModel& model);

    void implicitLoop
    (
        const scalar t0,
        const scalar dt,
        volScalarField& Vm,
        volScalarField& Iion,
        Field<Field<scalar>>& states,
        volScalarField& externalStimulusCurrent,
        const List<labelList>& stimulusCellIDsList,
        const List<scalar>& stimulusStartTimes,
        const scalar stimulusIntensity,
        const scalar stimulusDuration,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const volTensorField& conductivity,
        pimpleControl& pimple
    );
};

#endif
