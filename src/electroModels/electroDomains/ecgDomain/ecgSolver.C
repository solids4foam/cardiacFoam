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

#include "ecgSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam
{

defineTypeNameAndDebug(ECGSolver, 0);
defineRunTimeSelectionTable(ECGSolver, dictionary);


autoPtr<ECGSolver> ECGSolver::New(const dictionary& dict)
{
    const word solverType
    (
        dict.lookupOrDefault<word>("ecgSolver", "pseudoECG")
    );

    Info<< "Selecting ecgSolver " << solverType << nl;

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ecgSolver type " << solverType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ECGSolver>(ctorPtr(dict));
}

} // End namespace Foam

// ************************************************************************* //
