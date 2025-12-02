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
    tissue_(-1),
    solveVmWithinODESolver_(solveVmWithinODESolver)
{
    // Load exportedVariables list into variableExport_
    if (dict_.found("exportedVariables"))
    {
        dict_.lookup("exportedVariables") >> variableExport_;
    }
}

void Foam::ionicModel::setTissueFromDict()

{
     // --- First: check model override ---
     if (enforceModelTissue())
     {
         tissue_ = enforcedTissueFlag();
         Info<< "Tissue forced by model → " << enforcedTissueName() <<nl;
         return;
     }
    // --- Otherwise: read from dictionary ---
    bool hasTissue    = dict_.found("tissue");
    bool hasDimension = dict_.found("dimension");

    if (!hasTissue && !hasDimension)
    {
        FatalErrorInFunction
            << "Dictionary must specify either 'tissue' or 'dimension'."
            << exit(FatalError);
    }
    if (hasTissue && hasDimension)
    {
        FatalErrorInFunction
            << "Dictionary cannot specify BOTH 'tissue' and 'dimension'."
            << exit(FatalError);
    }
    word name;
    if (hasTissue)
    {
        dict_.lookup("tissue") >> name;

        if (!supportedTissueTypes().contains(name))
        {
            FatalErrorInFunction
                << "Unsupported tissue '" << name << "'. Allowed: "
                << supportedTissueTypes()
                << exit(FatalError);
        }

        tissue_ = (name=="epicardialCells")  ? 1 :
                  (name=="mCells")           ? 2 :
                  (name=="endocardialCells") ? 3 :
                  (name=="myocyte")          ? 4 : -1;
    }
    else // dimension
    {
        dict_.lookup("dimension") >> name;

        if (!supportedDimensions().contains(name))
        {
            FatalErrorInFunction
                << "Unsupported dimension '" << name << "'. Allowed: "
                << supportedDimensions()
                << exit(FatalError);
        }

        tissue_ = (name=="1D") ? 1 :
                  (name=="2D") ? 2 :
                  (name=="3D") ? 3 : -1;
    }

    Info<< "Selected " << name
        << " → internal tissue flag = " << tissue_ << nl;
}




// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ionicModel>
Foam::ionicModel::New
(
    const dictionary& dict,
    const label nIntegrationPoints,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
{
    const word modelType(dict.lookup("ionicModel"));
    

    Info<< "Selecting ionic model: " << modelType << endl;

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
        ctorPtr(dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
    );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModel::~ionicModel()
{}


// ************************************************************************* //