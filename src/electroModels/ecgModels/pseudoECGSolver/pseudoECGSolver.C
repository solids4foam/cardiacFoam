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

#include "pseudoECGSolver.H"

#include "ecgDomain.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(PseudoECGSolver, 0);
addToRunTimeSelectionTable(ECGSolver, PseudoECGSolver, dictionary);


void PseudoECGSolver::compute
(
    const ECGDomain& domain,
    scalarField& values
)
{
    // Gima-Rudy dipole:
    //   phi_pseudo(P) = -sum_c
    //     [ (conductivity . grad(Vm))_c . r_vec * V_c / |r|^3 ]
    //   r_vec = C_c - P

    const fvMesh& mesh = domain.mesh();
    const List<vector>& electrodePositions = domain.electrodePositions();

    const tmp<volVectorField> tgradVm = fvc::grad(domain.Vm());
    const vectorField& gradVm = tgradVm().primitiveField();

    const scalarField& volumes = mesh.V();
    const vectorField& cellCentres = mesh.C().primitiveField();
    const tensorField& conductivityField =
        domain.conductivity().primitiveField();

    const label nElectrodes = electrodePositions.size();

    values.setSize(nElectrodes);
    values = 0.0;

    forAll(cellCentres, cellI)
    {
        const vector dipole =
            (conductivityField[cellI] & gradVm[cellI]) * volumes[cellI];

        for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
        {
            const vector rVec =
                cellCentres[cellI] - electrodePositions[electrodeI];
            const scalar r = mag(rVec);

            if (r > VSMALL)
            {
                values[electrodeI] += (dipole & rVec)/(r*r*r);
            }
        }
    }

    for (label electrodeI = 0; electrodeI < nElectrodes; ++electrodeI)
    {
        reduce(values[electrodeI], sumOp<scalar>());
    }
}

} // End namespace Foam

// ************************************************************************* //
