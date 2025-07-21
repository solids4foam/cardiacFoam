/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ionicModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ionicModel, 0);
    defineRunTimeSelectionTable(ionicModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModel::ionicModel
(
    const dictionary& dict, const label num, const scalar initialDeltaT
)
:
    ODESystem(),
    odeSolver_(),
    dict_(dict),
    step_(num, initialDeltaT),
    tissue_(-1)  // initialize to invalid flag
{
    word tissueName;
    dict.lookup("tissue") >> tissueName;

    Info << "Tissue Name " << tissueName << endl;

    tissue_ = (tissueName == "epicardialCells") ? 1
            : (tissueName == "mCells")          ? 2
            : (tissueName == "endocardialCells") ? 3
            : -1;  // invalid flag

    if (tissue_ == -1)
    {
        FatalErrorInFunction
            << "Unknown tissue: " << tissueName << nl
            << exit(FatalError);
    }

    Info << "Tissue flag set to: " << tissue_ << endl;




}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ionicModel>
Foam::ionicModel::New
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const scalar initialDeltaT
)
{
    const word modelType(dict.lookup("ionicModel"));

    Info<< "Selecting ionic model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "ionicModel",
            modelType,
            *dictionaryConstructorTablePtr_
        )   << exit(FatalIOError);
    }

    return autoPtr<ionicModel>
    (
        ctorPtr(dict, nIntegrationPoints, initialDeltaT)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModel::~ionicModel()
{}


// ************************************************************************* //
