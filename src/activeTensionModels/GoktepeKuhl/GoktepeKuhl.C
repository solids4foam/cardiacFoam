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

#include "GoktepeKuhl.H"
#include "GoktepeKuhl_2004.H"
#include "activeTensionIO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GoktepeKuhl, 0);
    addToRunTimeSelectionTable(activeTensionModel, GoktepeKuhl, dictionary);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::GoktepeKuhl::GoktepeKuhl
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
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    currentVm_(0.0),
    currentAct_(0.0)
{
    Info<< nl << "Initialize GoktepeKuhl constants:" << nl;
    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI,     new scalarField(NUM_STATES, 0.0));
        ALGEBRAIC_.set(integrationPtI,  new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(integrationPtI,      new scalarField(NUM_STATES, 0.0));

        GoktepeKuhlinitConsts
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

void Foam::GoktepeKuhl::calculateTension
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
    const scalar tEnd   = (t + dt)  ;
    scalar step         = dt ;

    const CouplingSignalProvider& p = provider();
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& Ta_i  = STATESI[::Ta];
        scalar act    = p.signal(integrationPtI, CouplingSignal::Vm);
        scalar act2   = p.signal(integrationPtI, CouplingSignal::Act);
        currentVm_ = act;
        currentAct_ = act2;
        ALGEBRAICI[::AV_Vm] = act; // Convert to mV
        ALGEBRAICI[::AV_u] = act2;



        // Advance ODE system
        odeSolver_->solve(tStart, tEnd, STATESI, step);

        GoktepeKuhlcomputeVariables
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


void Foam::GoktepeKuhl::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);
    ALGEBRAIC_TMP[::AV_Vm] = currentVm_;
    ALGEBRAIC_TMP[::AV_u] = currentAct_;

    GoktepeKuhlcomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data()
    );

}


void Foam::GoktepeKuhl::jacobian
(
    const scalar,
    const scalarField&,
    scalarField&,
    scalarSquareMatrix&
) const
{
    notImplemented("Foam::GoktepeKuhl::jacobian(...)");
}

void Foam::GoktepeKuhl::debugPrintFields
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

    const label nAvailable = NUM_STATES + NUM_ALGEBRAIC + NUM_STATES + 1;
    wordList availableNames(nAvailable);
    scalarField values(nAvailable);

    for (label i = 0; i < NUM_STATES; ++i)
    {
        availableNames[i] = GoktepeKuhlSTATES_NAMES[i];
        values[i] = STATES_[cellI][i];
    }

    for (label i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        const label idx = NUM_STATES + i;
        availableNames[idx] = GoktepeKuhlALGEBRAIC_NAMES[i];
        values[idx] = ALGEBRAIC_[cellI][i];
    }

    for (label i = 0; i < NUM_STATES; ++i)
    {
        const label idx = NUM_STATES + NUM_ALGEBRAIC + i;
        availableNames[idx] =
            word("RATES_") + word(GoktepeKuhlSTATES_NAMES[i]);
        values[idx] = RATES_[cellI][i];
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

Foam::wordList Foam::GoktepeKuhl::availableFieldNames() const
{
    wordList names(NUM_STATES + NUM_ALGEBRAIC + NUM_STATES);

    for (label i = 0; i < NUM_STATES; ++i)
    {
        names[i] = GoktepeKuhlSTATES_NAMES[i];
    }

    for (label i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        names[NUM_STATES + i] = GoktepeKuhlALGEBRAIC_NAMES[i];
    }

    for (label i = 0; i < NUM_STATES; ++i)
    {
        names[NUM_STATES + NUM_ALGEBRAIC + i] =
            word("RATES_") + word(GoktepeKuhlSTATES_NAMES[i]);
    }

    return names;
}

void Foam::GoktepeKuhl::fieldValues(const label i, scalarField& values) const
{
    values.setSize(NUM_STATES + NUM_ALGEBRAIC + NUM_STATES);

    for (label s = 0; s < NUM_STATES; ++s)
    {
        values[s] = STATES_[i][s];
    }

    for (label a = 0; a < NUM_ALGEBRAIC; ++a)
    {
        values[NUM_STATES + a] = ALGEBRAIC_[i][a];
    }

    for (label r = 0; r < NUM_STATES; ++r)
    {
        values[NUM_STATES + NUM_ALGEBRAIC + r] = RATES_[i][r];
    }
}

// ************************************************************************* //
