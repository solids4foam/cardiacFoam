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
#include "IOdictionary.H"

#include "monodomainVerification/manufacturedFDAReference.H"

namespace Foam
{

//- Single source of truth for the manufactured displacement amplitude.
//  Both the boundary condition and the manufactured body force read it from
//  here so they can never disagree about which manufactured solution they are
//  imposing. The amplitude lives in the (top-level) electroMechanicalProperties
//  electromechanicalVerificationModel sub-dictionary.
inline vector manufacturedElectromechanicsAmplitude(const objectRegistry& db)
{
    IOdictionary emProps
    (
        IOobject
        (
            "electroMechanicalProperties",
            db.time().constant(),
            db.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word emModel(emProps.get<word>("electroMechanicalModel"));
    const dictionary& emCoeffs = emProps.subDict(emModel + "Coeffs");

    return emCoeffs.subDict("electromechanicalVerificationModel")
       .get<vector>("amplitude");
}


inline void computeManufacturedElectromechanicsVm
(
    scalarField& Vm,
    const vectorField& centres,
    const scalar t,
    const label dimension
)
{
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));
    
    computeManufacturedV(Vm, X, Y, Z, t, dimension);
}


inline void computeManufacturedElectromechanicsD
(
    vectorField& D,
    const vectorField& centres,
    const scalar t,
    const vector& amplitude
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
            amplitude.x()*x*x*y*st,
            amplitude.y()*y*y*z*st,
            amplitude.z()*z*z*x*st
        );
    }
}


inline scalar manufacturedElectromechanicsLambdaX
(
    const vector& X,
    const scalar t,
    const vector& amplitude
)
{
    const scalar st = Foam::sin(t);
    const scalar x = X.x();
    const scalar y = X.y();
    const scalar z = X.z();

    const scalar F11 = 1.0 + 2.0*amplitude.x()*x*y*st;
    const scalar F21 = 0.0;
    const scalar F31 = amplitude.z()*z*z*st;

    return Foam::sqrt(F11*F11 + F21*F21 + F31*F31);
}


inline void computeManufacturedElectromechanicsLambda
(
    scalarField& lambda,
    const vectorField& centres,
    const scalar t,
    const vector& amplitude
)
{
    lambda.setSize(centres.size());

    forAll(centres, cellI)
    {
        lambda[cellI] =
            manufacturedElectromechanicsLambdaX
            (
                centres[cellI],
                t,
                amplitude
            );
    }
}


inline void computeManufacturedElectromechanicsTa
(
    scalarField& Ta,
    const vectorField& centres,
    const scalar t,
    const vector& amplitude,
    const scalar Tmax,
    const scalar V0,
    const scalar gamma,
    const label dimension
)
{
    scalarField Vm;
    computeManufacturedElectromechanicsVm(Vm, centres, t, dimension);

    scalarField lambda;
    computeManufacturedElectromechanicsLambda
    (
        lambda,
        centres,
        t,
        amplitude
    );

    Ta.setSize(centres.size());

    const scalar V02 = V0*V0;

    forAll(centres, cellI)
    {
        const scalar Vm2 = Vm[cellI]*Vm[cellI];

        Ta[cellI] =
            Tmax
           *Vm2/(V02 + Vm2)
           *(1.0 + gamma*(lambda[cellI] - 1.0));
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
