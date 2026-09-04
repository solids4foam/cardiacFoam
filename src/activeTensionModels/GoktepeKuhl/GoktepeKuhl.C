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

void Foam::GoktepeKuhl::refreshRestartState(const fvMesh& mesh)
{
    forAll(STATES_, i)
    {
        scalarField& rates = RATES_[i];
        scalarField& algebraic = ALGEBRAIC_[i];
        const scalar drive = coupledDriveSignal(i);
        scalar u = (drive - CONSTANTS_[AC_Vr])/100.0;
        u = max(scalar(0.0), min(u, scalar(1.0)));
        algebraic[AV_Vm] = drive;
        algebraic[AV_u] = u;
        GoktepeKuhlcomputeVariables
        (
            mesh.time().value(),
            CONSTANTS_.data(),
            rates.data(),
            STATES_[i].data(),
            algebraic.data()
        );
    }
}

bool Foam::GoktepeKuhl::restartTension(scalarField& Ta) const
{
    if (Ta.size() != STATES_.size())
    {
        return false;
    }
    forAll(Ta, i)
    {
        Ta[i] = STATES_[i][::Ta];
    }
    return true;
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

    scalarField protoStates(NUM_STATES, 0.0);
    scalarField protoRates(NUM_STATES, 0.0);

    GoktepeKuhlinitConsts
    (
        CONSTANTS_.data(),
        protoRates.data(),
        protoStates.data()
    );

    if (dict.found("constants"))
    {
        const dictionary& cDict = dict.subDict("constants");
        for (label k = 0; k < NUM_CONSTANTS; ++k)
        {
            const word name(GoktepeKuhlCONSTANTS_NAMES[k]);
            if (cDict.found(name))
            {
                CONSTANTS_[k] = cDict.get<scalar>(name);
            }
        }
    }

    if (dict.found("initialStates"))
    {
        const dictionary& sDict = dict.subDict("initialStates");
        for (label k = 0; k < NUM_STATES; ++k)
        {
            const word name(GoktepeKuhlSTATES_NAMES[k]);
            if (sDict.found(name))
            {
                protoStates[k] = sDict.get<scalar>(name);
            }
        }
    }

    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI,    new scalarField(protoStates));
        ALGEBRAIC_.set(integrationPtI, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(integrationPtI,     new scalarField(NUM_STATES,    0.0));
    }
    Info<< CONSTANTS_ << nl;
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

    const scalar tStart = currentT_;
    const scalar tEnd   = currentT_ + currentDt_;
    scalar step         = currentDt_;

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
