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
#include "addToRunTimeSelectionTable.H"
#include "word.H"

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
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    currentDriveSignal_(0.0),
    driveSignal_(CouplingSignal::Act),
    driveSignalName_("Act")
{
    const word requestedSignal =
        dict_.lookupOrDefault<word>("couplingSignal", "Act");

    if (requestedSignal == "Act" || requestedSignal == "act")
    {
        driveSignal_     = CouplingSignal::Act;
        driveSignalName_ = "Act";
    }
    else if (requestedSignal == "Vm" || requestedSignal == "vm")
    {
        driveSignal_     = CouplingSignal::Vm;
        driveSignalName_ = "Vm";
    }
    else
    {
        FatalErrorInFunction
            << "Unknown GoktepeKuhl 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid options are: Act, Vm."
            << abort(FatalError);
    }

    Info<< nl << "Initialize GoktepeKuhl constants:" << nl;
    Info<< "GoktepeKuhl couplingSignal: " << driveSignalName_ << nl;

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

void Foam::GoktepeKuhl::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField&,
    scalarField& Ta
)
{
    if (Ta.size() != nIntegrationPoints_)
    {
        FatalErrorInFunction
            << "Ta.size() != nIntegrationPoints" << abort(FatalError);
    }

    const scalar tStart = t * 1000/12.9;
    const scalar tEnd   = (t + dt) * 1000/12.9;
    scalar step         = dt * 1000/12.9;

    const CouplingSignalProvider& p = provider();
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& Ta_i          = STATESI[::Ta];
        const scalar driveVal = p.signal(integrationPtI, driveSignal_);
        currentDriveSignal_   = driveVal;
        ALGEBRAICI[::AV_Vm]   = driveVal;
        ALGEBRAICI[::AV_u]    = driveVal;

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

        if (!debugPrintedNames().empty() && integrationPtI == monitorCell)
        {
            debugPrintFields(integrationPtI, tStart, tEnd, step);
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
    ALGEBRAIC_TMP[::AV_Vm] = currentDriveSignal_;
    ALGEBRAIC_TMP[::AV_u]  = currentDriveSignal_;

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
