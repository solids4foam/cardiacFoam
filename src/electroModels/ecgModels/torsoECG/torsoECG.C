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

#include "torsoECG.H"

#include "ecgDomain.H"
#include "ecgElectrodeCells.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(torsoECG, 0);
addToRunTimeSelectionTable(ecgSolver, torsoECG, dictionary);


torsoECG::torsoECG(const dictionary& dict)
:
    electrodeCells_(),
    cellsBuilt_(false)
{}


void torsoECG::buildElectrodeCells(const ecgDomain& domain)
{
    const volScalarField* phiEPtr = domain.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "torsoECG requires a global phiE field from the ECG "
            << "state provider."
            << exit(FatalError);
    }

    const fvMesh& phiEMesh = phiEPtr->mesh();
    const List<vector>& electrodes = domain.electrodePositions();

    electrodeCells_ = ecgElectrodeCells(phiEMesh, electrodes);

    Info<< "torsoECG: located " << electrodes.size()
        << " electrodes on phiE mesh '" << phiEMesh.name() << "'."
        << endl;

    cellsBuilt_ = true;
}


void torsoECG::solve
(
    ecgDomain& domain,
    scalar t0,
    scalar dt,
    scalarField& values
)
{
    (void)t0;
    (void)dt;

    if (!cellsBuilt_)
    {
        buildElectrodeCells(domain);
    }

    const volScalarField* phiEPtr = domain.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "torsoECG requires phiE, but domain.phiEPtr() returned null."
            << exit(FatalError);
    }

    const scalarField& phiE = phiEPtr->primitiveField();

    values.setSize(electrodeCells_.size());
    values = 0.0;

    forAll(electrodeCells_, electrodeI)
    {
        const label cellI = electrodeCells_[electrodeI];
        if (cellI >= 0)
        {
            values[electrodeI] = phiE[cellI];
        }
    }

    if (Pstream::parRun())
    {
        Pstream::listCombineGather(values, plusEqOp<scalar>());
        Pstream::broadcast(values);
    }
}

} // End namespace Foam

// ************************************************************************* //
