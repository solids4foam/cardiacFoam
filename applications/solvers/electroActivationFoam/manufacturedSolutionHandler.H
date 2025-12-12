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

#ifndef manufacturedSolutionHandler_H
#define manufacturedSolutionHandler_H

#include "ionicModel.H"


using namespace Foam;


class manufacturedSolutionHandler
{
private:
    const fvMesh& mesh_;
    ionicModel& ionicModel_;

    scalar dx_;
    int dim_;
    int N_;

public:
    // Constructor
    manufacturedSolutionHandler(const fvMesh& mesh, ionicModel& model);

    void initializeManufactured
    (
        volScalarField& Vm,
        PtrList<volScalarField>& outFields,
        scalar dxstructured,
        int dim
    );


    // Manufactured-solution error post-processing
    void postProcess
    (
        const volScalarField& Vm,
        const PtrList<volScalarField>& outFields,
        const scalar dt,
        const int nsteps,
        const bool solveExplicit
    ) const;


};

#endif
