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

#include "GoktepeKuhl.H"
#include "GoktepeKuhl_2004.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GoktepeKuhl, 0);
    addToRunTimeSelectionTable(activeTensionModel, GoktepeKuhl, dictionary);
}


// * * * * * * * * * * * * * * * * io* hooks  * * * * * * * * * * * * * * * //

const char* const* Foam::GoktepeKuhl::ioStateNames() const
{
    return GoktepeKuhlSTATES_NAMES;
}

const char* const* Foam::GoktepeKuhl::ioConstantNames() const
{
    return GoktepeKuhlCONSTANTS_NAMES;
}

const char* const* Foam::GoktepeKuhl::ioAlgebraicNames() const
{
    return GoktepeKuhlALGEBRAIC_NAMES;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    CONSTANTS_(NUM_CONSTANTS, 0.0)
{
    const word requestedSignal =
        dict_.lookupOrDefault<word>("couplingSignal", "Vm");

    if (!(requestedSignal == "Vm" || requestedSignal == "vm"))
    {
        FatalErrorInFunction
            << "Unknown GoktepeKuhl 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid option is: Vm."
            << abort(FatalError);
    }

    Info<< nl << "Initialize GoktepeKuhl constants:" << nl;
    Info<< "GoktepeKuhl couplingSignal: Vm" << nl;

    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI,    new scalarField(NUM_STATES,    0.0));
        ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(integrationPtI,     new scalarField(NUM_STATES,    0.0));

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

void Foam::GoktepeKuhl::solveAtPoint
(
    const label i,
    const scalar driveVal,
    const scalar /*lambda*/,
    scalar& Ta
)
{
    scalarField& STATESI    = STATES_[i];
    scalarField& ALGEBRAICI = ALGEBRAIC_[i];
    scalarField& RATESI     = RATES_[i];

    scalar& Ta_i = STATESI[::Ta];

    // Map Vm (mV) to Aliev-Panfilov activation variable u in [0,1].
    scalar uSignal = (driveVal - CONSTANTS_[AC_Vr]) / 100.0;

    uSignal = max(scalar(0.0), min(uSignal, scalar(1.0)));

    ALGEBRAICI[::AV_Vm] = driveVal;
    ALGEBRAICI[::AV_u]  = uSignal;

    const scalar tStart = currentT_ * 1000/100;
    const scalar tEnd   = (currentT_ + currentDt_) * 1000/100;
    scalar step         = currentDt_ * 1000/100;

    odeSolver_->solve(tStart, tEnd, STATESI, step);

    GoktepeKuhlcomputeVariables
    (
        tEnd,
        CONSTANTS_.data(),
        RATESI.data(),
        STATESI.data(),
        ALGEBRAICI.data()
    );

    Ta = Ta_i;

    if (!debugPrintedNames().empty() && i == 0)
    {
        debugPrintFields(i, tStart, tEnd, step);
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
    scalar uSignal = (currentDriveSignal_ - CONSTANTS_[AC_Vr]) / 100.0;

    uSignal = max(scalar(0.0), min(uSignal, scalar(1.0)));

    ALGEBRAIC_TMP[::AV_Vm] = currentDriveSignal_;
    ALGEBRAIC_TMP[::AV_u]  = uSignal;

    GoktepeKuhlcomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
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

// ************************************************************************* //
