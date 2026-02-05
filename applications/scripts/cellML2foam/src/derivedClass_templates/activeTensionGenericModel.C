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
    STATES_(num),
    ALGEBRAIC_(num),
    RATES_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0)
{
    Info<< nl << "Initialize activeTensionGenericModel constants:" << nl;
    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI,    new scalarField(NUM_STATES, 0.0));
        ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(integrationPtI,     new scalarField(NUM_STATES, 0.0));

        activeTensionGenericModelinitConsts
        (
            CONSTANTS_.data(),
            RATES_[integrationPtI].data(),
            STATES_[integrationPtI].data()
        );
    }
    Info<< CONSTANTS_ << nl;

    label i0 = rand() % STATES_.size();
    Info<< "initial states:" << nl;
    Info<< STATES_[i0] << nl;
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

    const CouplingSignalProvider& p = provider();
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& Ta_i  = STATESI[::Ta];
        scalar& act   = STATESI[::u];

        act = p.signal(integrationPtI, CouplingSignal::Act);

        // Advance ODE system
        odeSolver_->solve(tStart, tEnd, STATESI, step);

        activeTensionGenericModelcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data()
        );

        Ta[integrationPtI] = Ta_i;

        const wordList printedNames = debugPrintedNames();
        if (!printedNames.empty() && integrationPtI == monitorCell)
        {
            debugPrintFields(integrationPtI, tStart, tEnd, step, act);
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
    scalarField rates(NUM_STATES, 0.0);

    activeTensionGenericModelcomputeVariables
    (
        t,
        CONSTANTS_.data(),
        rates.data(),
        const_cast<scalarField&>(y).data(),
        alg.data()
    );

    dydt[::Ta] = rates[::Ta];
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

void Foam::activeTensionGenericModel::debugPrintFields
(
    const label cellI,
    const scalar t1,
    const scalar t2,
    const scalar step,
    const scalar act
) const
{
    const wordList printedNames = debugPrintedNames();
    if (printedNames.empty())
    {
        return;
    }

    const label nAvailable = NUM_STATES + NUM_ALGEBRAIC + 1;
    wordList availableNames(nAvailable);
    scalarField values(nAvailable);

    for (label i = 0; i < NUM_STATES; ++i)
    {
        availableNames[i] = activeTensionGenericModelSTATES_NAMES[i];
        values[i] = STATES_[cellI][i];
    }

    for (label i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        const label idx = NUM_STATES + i;
        availableNames[idx] = activeTensionGenericModelALGEBRAIC_NAMES[i];
        values[idx] = ALGEBRAIC_[cellI][i];
    }

    availableNames[nAvailable - 1] = "Act";
    values[nAvailable - 1] = act;

    activeTensionIO::debugPrintFields
    (
        printedNames,
        availableNames,
        values,
        cellI,
        t1,
        t2,
        step
    );
}

Foam::wordList Foam::activeTensionGenericModel::availableFieldNames() const
{
    wordList names(NUM_STATES + NUM_ALGEBRAIC);

    for (label i = 0; i < NUM_STATES; ++i)
    {
        names[i] = activeTensionGenericModelSTATES_NAMES[i];
    }

    for (label i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        names[NUM_STATES + i] = activeTensionGenericModelALGEBRAIC_NAMES[i];
    }

    return names;
}

void Foam::activeTensionGenericModel::fieldValues
(
    const label i,
    scalarField& values
) const
{
    values.setSize(NUM_STATES + NUM_ALGEBRAIC);

    for (label s = 0; s < NUM_STATES; ++s)
    {
        values[s] = STATES_[i][s];
    }

    for (label a = 0; a < NUM_ALGEBRAIC; ++a)
    {
        values[NUM_STATES + a] = ALGEBRAIC_[i][a];
    }
}

// ************************************************************************* //
