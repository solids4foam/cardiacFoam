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

#include "ionicModelFDA.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ionicModelFDA, 0);
    defineRunTimeSelectionTable(ionicModelFDA, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModelFDA::ionicModelFDA
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

void Foam::ionicModelFDA::setTissueFromDict()
{
    word tissueName;
    dict_.lookup("dimension") >> tissueName;
    

    List <word> tissues = supportedTissues();


    tissue_ = (tissueName == "1D") ? 1
            : (tissueName == "2D") ? 2
            : (tissueName == "3D") ? 3
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

Foam::autoPtr<Foam::ionicModelFDA>
Foam::ionicModelFDA::New
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

    return autoPtr<ionicModelFDA>
    (
        ctorPtr(dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModelFDA::~ionicModelFDA()
{}


// ************************************************************************* //
