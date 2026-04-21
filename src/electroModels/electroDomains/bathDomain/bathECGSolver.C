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

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  Stub for future bath-ECG coupling solver.
\*---------------------------------------------------------------------------*/

#include "bathECGSolver.H"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"

namespace Foam
{

defineTypeNameAndDebug(BathECGSolver, 0);
defineRunTimeSelectionTable(BathECGSolver, dictionary);


autoPtr<BathECGSolver> BathECGSolver::New(const dictionary& dict)
{
    const word solverType
    (
        dict.lookupOrDefault<word>("bathSolver", "bidomainBathECG")
    );

    Info<< "Selecting bathECGSolver " << solverType << nl;

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown bathECGSolver type " << solverType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<BathECGSolver>(ctorPtr(dict));
}

} // End namespace Foam

// ************************************************************************* //
