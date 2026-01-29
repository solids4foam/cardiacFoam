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

#include <cmath>

#include "activeTensionGenericModel.H"
#include "activeTensionGenericModel_YYYY.H"
#include "activeTensionIO.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(activeTensionGenericModel, 0);
    addToRunTimeSelectionTable(activeTensionModel, activeTensionGenericModel, dictionary);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::activeTensionGenericModel::activeTensionGenericModel
(
    const dictionary& dict,
    const label num
)
:
    activeTensionModel(dict, num),
    ODESystem(),
    odeSolver_(ODESolver::New(*this, dict_)),
    internalVariables_(num),
    algebraicVariables_(num),
    constants_(NUM_CONSTANTS, 0.0)
{
    forAll(internalVariables_, integrationPtI)
    {
        internalVariables_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));
        algebraicVariables_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));

        scalarField rates(NUM_STATES, 0.0);
        activeTensionGenericModelinitConsts
        (
            constants_.data(),
            rates.data(),
            internalVariables_[integrationPtI].data()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activeTensionGenericModel::calculateTension
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

    const scalar tStart = t;
    const scalar tEnd   = t + dt;
    scalar step         = dt;

    const label monitorCell = 0;

    forAll(internalVariables_, integrationPtI)
    {
        scalarField& yStart = internalVariables_[integrationPtI];
        scalarField& alg    = algebraicVariables_[integrationPtI];

        // Advance ODE system
        odeSolver_->solve(tStart, tEnd, yStart, step);

        scalarField rates(NUM_STATES, 0.0);
        activeTensionGenericModelcomputeVariables
        (
            tEnd,
            constants_.data(),
            rates.data(),
            yStart.data(),
            alg.data(),
            false
        );

        Ta[integrationPtI] = alg[AV_Tension];

        const wordList printedNames = debugPrintedNames();
        if (!printedNames.empty() && integrationPtI == monitorCell)
        {
            wordList availableNames(2);
            availableNames[0] = "Ta";
            availableNames[1] = "Act";

            scalarField values(2);
            values[0] = Ta[integrationPtI];
            values[1] = Act[integrationPtI];

            activeTensionIO::debugPrintFields
            (
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


void Foam::activeTensionGenericModel::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField alg(NUM_ALGEBRAIC, 0.0);

    activeTensionGenericModelcomputeVariables
    (
        t,
        constants_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        alg.data(),
        false
    );
}


void Foam::activeTensionGenericModel::jacobian
(
    const scalar,
    const scalarField&,
    scalarField&,
    scalarSquareMatrix&
) const
{
    notImplemented("Foam::activeTensionGenericModel::jacobian(...)");
}

// ************************************************************************* //
