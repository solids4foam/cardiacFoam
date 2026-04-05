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

#include "explicit1DSolver.H"
#include "graphConductionSystemDomain.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(Explicit1DSolver, 0);
addToRunTimeSelectionTable
(
    GraphConductionSystemSolver,
    Explicit1DSolver,
    dictionary
);


void Explicit1DSolver::advance
(
    GraphConductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    scalarField& Vm = domain.membranePotential();
    scalarField& Iion = domain.ionicCurrent();

    domain.ionicModelRef().solveODE(t0, dt, Vm, Iion);

    laplacianBuffer_.setSize(Vm.size());
    laplacianBuffer_ = 0.0;

    const labelList& edgeA = domain.edgeStartNodes();
    const labelList& edgeB = domain.edgeEndNodes();
    const scalarField& edgeLength = domain.edgeLengths();
    const scalarField& edgeConductance = domain.edgeConductances();

    forAll(edgeA, edgeI)
    {
        const label A = edgeA[edgeI];
        const label B = edgeB[edgeI];
        const scalar flux =
            edgeConductance[edgeI]
          * (Vm[B] - Vm[A])
          / sqr(edgeLength[edgeI]);

        laplacianBuffer_[A] += flux;
        laplacianBuffer_[B] -= flux;
    }

    appliedCurrentBuffer_.setSize(Vm.size());
    domain.assembleAppliedCurrent(t0, appliedCurrentBuffer_);

    const scalar dtOverChiCm = dt / (domain.chi() * domain.Cm());
    forAll(Vm, nodeI)
    {
        Vm[nodeI] += dtOverChiCm
                   * (laplacianBuffer_[nodeI] + appliedCurrentBuffer_[nodeI])
                   - dt * Iion[nodeI];
    }

    domain.reportAdvanceDiagnostics(t0, dt);
}

} // End namespace Foam

// ************************************************************************* //
