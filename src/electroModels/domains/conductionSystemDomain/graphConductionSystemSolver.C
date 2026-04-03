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

#include "graphConductionSystemSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam
{

defineTypeNameAndDebug(GraphConductionSystemSolver, 0);
defineRunTimeSelectionTable(GraphConductionSystemSolver, dictionary);


autoPtr<GraphConductionSystemSolver> GraphConductionSystemSolver::New
(
    const dictionary& dict
)
{
    const word solverType
    (
        dict.lookupOrDefault<word>
        (
            "ConductionSystemSolver",
            "MonoDomain1DSolver"
        )
    );

    Info<< "Selecting ConductionSystemSolver " << solverType << nl;

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ConductionSystemSolver type " << solverType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GraphConductionSystemSolver>(ctorPtr(dict));
}

} // End namespace Foam

// ************************************************************************* //
