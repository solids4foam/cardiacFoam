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

#include "newionicModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newionicModel, 0);
    defineRunTimeSelectionTable(newionicModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newionicModel::newionicModel
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ODESystem(),
    odeSolver_(),
    dict_(dict),
    step_(num, initialDeltaT),
    tissue_(-1),  // initialize to invalid flag
    solveVmWithinODESolver_(solveVmWithinODESolver)
{
    setTissueFromDict(); 

}

void Foam::newionicModel::setTissueFromDict()
{
    word tissueName;
    dict_.lookup("tissue") >> tissueName;
    

    List <word> tissues = supportedTissues();

    if (tissues.size() == 1)
    {
        tissueName = tissues[0];
        Info << "Atrial model or ventricular myocyte model with no endo-mcell-epi distinction "
             << tissueName << ". Overriding input tissue." << endl;
    }

    tissue_ = (tissueName == "epicardialCells") ? 1
            : (tissueName == "mCells")          ? 2
            : (tissueName == "endocardialCells") ? 3
            : (tissueName == "myocyte") ? 4
            : -1;  // invalid flag

    if (tissue_ == -1 || !tissues.contains(tissueName))
    {
        FatalErrorInFunction
            << "Unknown or unsupported tissue: " << tissueName << nl
            << "Supported tissues are: " << tissues << nl
            << exit(FatalError);
    }
    Info << "Tissue Name set to: " << tissueName << endl;
    Info << "Tissue flag set to: " << tissue_ << endl;
    
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::newionicModel>
Foam::newionicModel::New
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
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

    return autoPtr<newionicModel>
    (
        ctorPtr(dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newionicModel::~newionicModel()
{}


// ************************************************************************* //