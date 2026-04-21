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
  Stub for future bidomain-bath ECG coupling solver.
\*---------------------------------------------------------------------------*/

#include "bidomainBathECGSolver.H"

#include "bathDomain.H"
#include "fvScalarMatrix.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

addToRunTimeSelectionTable(BathECGSolver, bidomainBathECGSolver, dictionary);


bidomainBathECGSolver::bidomainBathECGSolver(const dictionary& dict)
:
    phiEReferenceCell_
    (
        dict.lookupOrDefault<label>("phiEReferenceCell", 0)
    ),
    phiEReferenceValue_
    (
        dict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0)
    ),
    reportBathExtrema_
    (
        dict.lookupOrDefault<Switch>("reportBathExtrema", false)
    )
{}


void bidomainBathECGSolver::solve
(
    BathDomain& domain,
    scalar t0,
    scalar dt
)
{
    volScalarField& phiE = domain.phiERef();
    const volScalarField& sigmaBath = domain.sigmaBath();

    fvScalarMatrix phiEqn
    (
        fvm::laplacian(sigmaBath, phiE)
     == -domain.interfaceSourceField()
    );

    phiEqn.setReference(phiEReferenceCell_, phiEReferenceValue_);

    phiEqn.solve();

    if (reportBathExtrema_)
    {
        Info<< "Bath phiE: min = " << min(phiE).value()
            << "  max = " << max(phiE).value() << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
