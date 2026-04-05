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
    const word& solverType,
    const dictionary& coeffs
)
{
    Info<< nl << "Selecting myocardiumSolver " << solverType << endl;

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            coeffs,
            "myocardiumSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<myocardiumSolver>(ctorPtr(mesh, coeffs));
}

} // End namespace Foam

// ************************************************************************* //
