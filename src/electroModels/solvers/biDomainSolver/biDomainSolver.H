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
    BiDomainSolver

Description
    Generic solver for cardiac electrophysiology based on the bidomain
    reaction-diffusion system.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    BiDomainSolver.C
\*---------------------------------------------------------------------------*/

#ifndef BiDomainSolver_H
#define BiDomainSolver_H

#include "electroModel.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{
namespace electroModels
{

class BiDomainSolver
:
    public electroModel
{
    volScalarField phiE_;

    BiDomainSolver(const BiDomainSolver&);
    void operator=(const BiDomainSolver&);

protected:

    autoPtr<ionicModel> ionicModelPtr_;
    autoPtr<electroVerificationModel> verificationModelPtr_;
    PtrList<volScalarField> outFields_;
    wordList postProcessFieldNames_;
    PtrList<volScalarField> postProcessFields_;

    void configureSystem();

    BiDomainSolver
    (
        const word& type,
        Time& runTime,
        const word& region = dynamicFvMesh::defaultRegion
    );

public:
    OverrideTypeName("BiDomainSolver");

    BiDomainSolver
    (
        Time& runTime,
        const word& region = dynamicFvMesh::defaultRegion
    );

    virtual ~BiDomainSolver() {}
};

} // End namespace electroModels
} // End namespace Foam

#endif

// ************************************************************************* //
