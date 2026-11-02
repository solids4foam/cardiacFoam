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
#include "torsoCoupling.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(torsoECG, 0);
addToRunTimeSelectionTable(ecgSolver, torsoECG, dictionary);


torsoECG::torsoECG(const dictionary& dict)
:
    leadFieldPtr_(),
    electrodeCells_(),
    cellsBuilt_(false)
{
    const dictionary& electroProperties = dict.topDict();

    if
    (
        torsoCoupling::read(electroProperties)
     == torsoCoupling::type::leadField
    )
    {
        leadFieldPtr_.reset
        (
            new leadFieldECG(electroProperties, dict.dictName())
        );
    }
}


void torsoECG::sample(const ecgDomain& domain, scalarField& values)
{
    const volScalarField& phiEField = domain.phiE();

    if (!cellsBuilt_)
    {
        electrodeCells_ = ecgElectrodeCells
        (
            phiEField.mesh(),
            domain.electrodeNames(),
            domain.electrodePositions()
        );
        cellsBuilt_ = true;

        Info<< "torsoECG: located " << electrodeCells_.size()
            << " electrodes on phiE mesh '" << phiEField.mesh().name()
            << "'." << endl;
    }

    const scalarField& phiE = phiEField.primitiveField();

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

    if (!leadFieldPtr_.valid())
    {
        sample(domain, values);
        return;
    }

    if (!leadFieldPtr_->loaded())
    {
        leadFieldPtr_->load
        (
            domain, domain.electrodeNames(), domain.electrodePositions()
        );
    }

    leadFieldPtr_->evaluate(domain.Vm(), values);
}

} // End namespace Foam

// ************************************************************************* //
