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


#include "NashPanfilov.H"
#include "activeTensionIO.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(NashPanfilov, 0);
    addToRunTimeSelectionTable(activeTensionModel, NashPanfilov, dictionary);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::NashPanfilov::NashPanfilov
(
    const dictionary& dict,
    const label num
)
:
    activeTensionModel(dict, num),
    ODESystem(),
    odeSolver_(ODESolver::New(*this, dict_)),
    internalVariables_(num),
    kTa_(dict.lookupOrDefault<scalar>("kTa", 47.9))
{
    // Lookup initial Ta
    const scalar initialTa(dict.lookupOrDefault<scalar>("initialTa", 0.0143524));

    forAll(internalVariables_, integrationPtI)
    {
        // y = [Ta, act]
        internalVariables_.set(integrationPtI, new scalarField(nEqns() + 1, 0.0));
        internalVariables_[integrationPtI][0] = initialTa; // Ta
        internalVariables_[integrationPtI][1] = 0.0;       // act
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::NashPanfilov::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField& Act,
    scalarField& Ta
)
{
    if (Ta.size() != nIntegrationPoints_)
    {
        FatalErrorInFunction
            << "Ta.size() != nIntegrationPoints" << abort(FatalError);
    }

    // (Optional) ensure provider exists + required signals exist
    // Usually call this once from the driver, but safe to call here too.
    // validateProvider();

    // interpret t as "current time", integrate from t-dt to t
    const scalar tStart = t;
    const scalar tEnd   = t + dt;
    scalar step         = dt ;

    const CouplingSignalProvider& p = provider();
    const label monitorCell = 0;

    forAll(internalVariables_, integrationPtI)
    {
        scalarField& yStart = internalVariables_[integrationPtI];

        // y = [Ta, act]
        scalar& Ta_i  = yStart[0];
        scalar& act   = yStart[1];

        act = p.signal(integrationPtI, CouplingSignal::Act);

        // Advance ODE system (Ta evolves, act is held constant via dydt[1]=0)
        odeSolver_->solve(tStart, tEnd, yStart, step);

        Ta[integrationPtI] = Ta_i;

        const wordList printedNames = debugPrintedNames();
        if (!printedNames.empty() && integrationPtI == monitorCell)
        {
            wordList availableNames(2);
            availableNames[0] = "Ta";
            availableNames[1] = "Act";
            scalarField values(2);
            values[0] = Ta_i;
            values[1] = act;
            activeTensionIO::debugPrintFields(
                printedNames,
                availableNames,
                values,
                integrationPtI,
                tStart,
                tEnd,
                step
            );
        }
    }
}


void Foam::NashPanfilov::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    const scalar& Ta  = y[0];
    const scalar& act = y[1];

    scalar& dTadt = dydt[0];

    //What do I do to the state rate that comes from ionicMOdel ??? ls
    dydt[1] = 0.0;

    scalar epsilonU = (act > 0.05 ? 10.0 : 1.0);

    dTadt = epsilonU*(kTa_*act - Ta);

}

void Foam::NashPanfilov::jacobian
(
    const scalar,
    const scalarField&,
    scalarField&,
    scalarSquareMatrix&
) const
{
    notImplemented("Foam::NashPanfilov::jacobian(...)");
}
