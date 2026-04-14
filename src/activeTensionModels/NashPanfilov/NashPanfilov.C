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
#include "NashPanfilov_2004.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NashPanfilov, 0);
    addToRunTimeSelectionTable(activeTensionModel, NashPanfilov, dictionary);
}


// * * * * * * * * * * * * * * * * io* hooks  * * * * * * * * * * * * * * * //

const char* const* Foam::NashPanfilov::ioStateNames() const
{
    return NashPanfilovSTATES_NAMES;
}

const char* const* Foam::NashPanfilov::ioConstantNames() const
{
    return NashPanfilovCONSTANTS_NAMES;
}

const char* const* Foam::NashPanfilov::ioAlgebraicNames() const
{
    return NashPanfilovALGEBRAIC_NAMES;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NashPanfilov::NashPanfilov
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
            << "Unknown NashPanfilov 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid option is: Vm."
            << abort(FatalError);
    }

    Info<< nl << "Initialize NashPanfilov constants:" << nl;
    Info<< "NashPanfilov couplingSignal: Vm" << nl;
    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI,    new scalarField(NUM_STATES,    0.0));
        ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(integrationPtI,     new scalarField(NUM_STATES,    0.0));

        NashPanfilovinitConsts
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

void Foam::NashPanfilov::solveAtPoint
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

    // Map Vm (mV) to normalized activation variable u in [0,1].
    scalar uSignal =
        (driveVal - CONSTANTS_[AC_Vr])
      / (CONSTANTS_[AC_Vp] - CONSTANTS_[AC_Vr]);

    uSignal = max(scalar(0.0), min(uSignal, scalar(1.0)));

    ALGEBRAICI[::AV_u] = uSignal;

    // Nash-Panfilov ODE was derived in the same normalized time used by
    // Aliev-Panfilov (t* = t_ms / 12.9), so convert OpenFOAM seconds first.
    const scalar tStart = currentT_ * 1000.0 / 12.9;
    const scalar tEnd   = (currentT_ + currentDt_) * 1000.0 / 12.9;
    scalar step         = currentDt_ * 1000.0 / 12.9;

    odeSolver_->solve(tStart, tEnd, STATESI, step);

    NashPanfilovcomputeVariables
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


void Foam::NashPanfilov::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);
    scalar uSignal =
        (currentDriveSignal_ - CONSTANTS_[AC_Vr])
      / (CONSTANTS_[AC_Vp] - CONSTANTS_[AC_Vr]);

    uSignal = max(scalar(0.0), min(uSignal, scalar(1.0)));
    ALGEBRAIC_TMP[::AV_u] = uSignal;

    NashPanfilovcomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALGEBRAIC_TMP.data()
    );
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

// ************************************************************************* //
