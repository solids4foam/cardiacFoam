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
    manufacturedElectromechanicsVerifier

Description
    Manufactured electromechanics verifier.

SourceFiles
    manufacturedElectromechanicsVerifier.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedElectromechanicsVerifier_H
#define manufacturedElectromechanicsVerifier_H

#include "electromechanicalVerificationModel.H"
#include "overrideTypeName.H"
#include "Switch.H"

namespace Foam
{

class manufacturedElectromechanicsVerifier
:
    public electromechanicalVerificationModel
{
    scalar uAmplitude_;
    scalar vAmplitude_;
    scalar wAmplitude_;
    scalar Tmax_;
    scalar V0_;
    scalar gamma_;
    scalar TaScale_;
    Switch initializeFields_;
    Switch enforceExactFields_;
    mutable bool errorsReported_;

    void setExactFields(volScalarField& Vm, volVectorField& D) const;

public:

    OverrideTypeName("manufacturedElectromechanicsVerifier");

    explicit manufacturedElectromechanicsVerifier(const dictionary& dict);

    virtual void initialize
    (
        volScalarField& Vm,
        volVectorField& D
    ) override;

    virtual void preSolve
    (
        volScalarField& Vm,
        volVectorField& D
    ) override;

    virtual bool shouldPostProcess
    (
        const volScalarField& Vm,
        const volVectorField& D
    ) const override;

    virtual void postProcess
    (
        const volScalarField& Vm,
        const volVectorField& D,
        const volScalarField& Ta
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
