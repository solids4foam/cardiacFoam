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

#include "myocardiumSolver.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(myocardiumSolver, 0);
defineRunTimeSelectionTable(myocardiumSolver, dictionary);


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

autoPtr<myocardiumSolver> myocardiumSolver::New
(
    const fvMesh& mesh,
    const dictionary& electroProperties
)
{
    const word modelType(electroProperties.lookup("myocardiumSolver"));

    Info<< nl << "Selecting myocardiumSolver " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            electroProperties,
            "myocardiumSolver",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("myocardiumSolver::New")
            << "Unknown myocardiumSolver type " << modelType
            << nl << nl
            << "Valid types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<myocardiumSolver>(ctorPtr(mesh, electroProperties));
}


autoPtr<myocardiumSolver> myocardiumSolver::NewBidomain
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    volScalarField& phiE
)
{
    // The bidomain solver type name is still looked up from the dict for
    // consistency, but we pass phiE through a separate construction path.
    // BiDomainDiffusionSolver registers itself under "bidomainSolver" via
    // addToRunTimeSelectionTable, but its constructor takes (mesh, dict, phiE).
    // We create it directly here to avoid duplicating the factory machinery.
    // Other solver types requested here are a FatalError.

    const word modelType(electroProperties.lookup("myocardiumSolver"));

    if (modelType != "bidomainSolver")
    {
        FatalErrorInFunction
            << "NewBidomain() called but myocardiumSolver is '"
            << modelType << "', not 'bidomainSolver'."
            << exit(FatalError);
    }

    // BiDomainDiffusionSolver is forward-declared here; the .C file that
    // includes it handles the actual construction.
    // Return via the standard table to keep the factory pattern clean.
    return New(mesh, electroProperties);
}

} // End namespace Foam

// ************************************************************************* //
