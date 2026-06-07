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
    manufacturedElectromechanicsReference

Description
    Manufactured electromechanics reference fields.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedElectromechanicsReference_H
#define manufacturedElectromechanicsReference_H

#include "mathematicalConstants.H"
#include "scalarField.H"
#include "vectorField.H"

namespace Foam
{

inline void computeManufacturedElectromechanicsVm
(
    scalarField& Vm,
    const vectorField& centres,
    const scalar t
)
{
    Vm.setSize(centres.size());

    const scalar pi = constant::mathematical::pi;
    const scalar ct = Foam::cos(t);

    forAll(centres, cellI)
    {
        const scalar x = centres[cellI].x();
        const scalar y = centres[cellI].y();
        const scalar z = centres[cellI].z();

        Vm[cellI] =
            Foam::sin(pi*x)
           *Foam::sin(pi*y)
           *Foam::sin(pi*z)
           *ct;
    }
}


inline void computeManufacturedElectromechanicsD
(
    vectorField& D,
    const vectorField& centres,
    const scalar t,
    const scalar uAmplitude,
    const scalar vAmplitude,
    const scalar wAmplitude
)
{
    D.setSize(centres.size());

    const scalar st = Foam::sin(t);

    forAll(centres, cellI)
    {
        const scalar x = centres[cellI].x();
        const scalar y = centres[cellI].y();
        const scalar z = centres[cellI].z();

        D[cellI] = vector
        (
            uAmplitude*x*x*y*st,
            vAmplitude*y*y*z*st,
            wAmplitude*z*z*x*st
        );
    }
}


inline scalar manufacturedElectromechanicsLambdaX
(
    const vector& X,
    const scalar t,
    const scalar uAmplitude,
    const scalar wAmplitude
)
{
    const scalar st = Foam::sin(t);
    const scalar x = X.x();
    const scalar y = X.y();
    const scalar z = X.z();

    const scalar F11 = 1.0 + 2.0*uAmplitude*x*y*st;
    const scalar F21 = 0.0;
    const scalar F31 = wAmplitude*z*z*st;

    return Foam::sqrt(F11*F11 + F21*F21 + F31*F31);
}


inline void computeManufacturedElectromechanicsTa
(
    scalarField& Ta,
    const vectorField& centres,
    const scalar t,
    const scalar uAmplitude,
    const scalar wAmplitude,
    const scalar Tmax,
    const scalar V0,
    const scalar gamma
)
{
    scalarField Vm;
    computeManufacturedElectromechanicsVm(Vm, centres, t);

    Ta.setSize(centres.size());

    const scalar V02 = V0*V0;

    forAll(centres, cellI)
    {
        const scalar Vm2 = Vm[cellI]*Vm[cellI];
        const scalar lambda =
            manufacturedElectromechanicsLambdaX
            (
                centres[cellI],
                t,
                uAmplitude,
                wAmplitude
            );

        Ta[cellI] =
            Tmax
           *Vm2/(V02 + Vm2)
           *(1.0 + gamma*(lambda - 1.0));
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
