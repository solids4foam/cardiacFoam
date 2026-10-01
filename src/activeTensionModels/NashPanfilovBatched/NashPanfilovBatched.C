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

#include "NashPanfilovBatched.H"
#include "NashPanfilovBatch.H"
#include "addToRunTimeSelectionTable.H"
#include "../NashPanfilov/NashPanfilov_2004Names.H"

#ifdef HAS_CUDA
namespace Foam
{
void launchNashPanfilovBatchKernel
(
    const double* driveSignals,
    const double* constants,
    int nCells,
    const double* states,
    double* rates,
    double* algebraics
);
}
#endif

namespace Foam
{
    defineTypeNameAndDebug(NashPanfilovBatched, 0);
    addToRunTimeSelectionTable(activeTensionModel, NashPanfilovBatched, dictionary);
}


// * * * * * * * * * * * * * * * * io* hooks  * * * * * * * * * * * * * * * //

const char* const* Foam::NashPanfilovBatched::ioStateNames() const
{
    static const char* const names[NUM_STATES] = { "Ta" };
    return names;
}

const char* const* Foam::NashPanfilovBatched::ioConstantNames() const
{
    static const char* const names[NUM_CONSTANTS] = {
        "AC_Vp", "AC_Vr", "AC_Vth", "AC_e0", "AC_kTa"
    };
    return names;
}

const char* const* Foam::NashPanfilovBatched::ioAlgebraicNames() const
{
    static const char* const names[NUM_ALGEBRAIC] = {
        "AV_u", "AV_e"
    };
    return names;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NashPanfilovBatched::NashPanfilovBatched
(
    const dictionary& dict,
    const label num
)
:
    batchedActiveTensionModel(dict, num, NUM_STATES, NUM_ALGEBRAIC),
    CONSTANTS_(NUM_CONSTANTS, 0.0)
{
    const word requestedSignal = dict_.lookupOrDefault<word>("couplingSignal", "Vm");

    if (!(requestedSignal == "Vm" || requestedSignal == "vm"))
    {
        FatalErrorInFunction
            << "Unknown NashPanfilovBatched 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid option is: Vm."
            << abort(FatalError);
    }

    Info<< nl << "Initialize NashPanfilovBatched constants:" << nl;

    // Default Initialization
    CONSTANTS_[AC_Vp] = 20.0;
    CONSTANTS_[AC_Vr] = -80.0;
    CONSTANTS_[AC_Vth] = -70.0;
    CONSTANTS_[AC_e0] = 1.0;
    CONSTANTS_[AC_kTa] = 47.9;

    scalarField protoStates(NUM_STATES, 0.0);

    if (dict.found("constants"))
    {
        const dictionary& cDict = dict.subDict("constants");
        for (label k = 0; k < NUM_CONSTANTS; ++k)
        {
            const word name(ioConstantNames()[k]);
            if (cDict.found(name)) CONSTANTS_[k] = cDict.get<scalar>(name);
        }
    }

    if (dict.found("initialStates"))
    {
        const dictionary& sDict = dict.subDict("initialStates");
        for (label k = 0; k < NUM_STATES; ++k)
        {
            const word name(ioStateNames()[k]);
            if (sDict.found(name)) protoStates[k] = sDict.get<scalar>(name);
        }
    }

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = protoStates[stateI];
        }
    }

    Info<< CONSTANTS_ << nl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::NashPanfilovBatched::evaluateHotPathStateForCell
(
    const label cellI,
    const scalar modelTime,
    const scalar driveSignal,
    const scalar lambda,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    NashPanfilovComputeVariablesBatch
    (
        driveSignal,
        CONSTANTS_.data(),
        1,
        0,
        const_cast<scalar*>(stateValues.begin()),
        rateValues.begin(),
        algebraicValues.begin()
    );
}


#ifdef HAS_CUDA
void Foam::NashPanfilovBatched::prepareCUDAInputs
(
    const scalarField& driveSignals,
    const scalarField& lambda
) const
{
    // NashPanfilov receives Vm through its dedicated device drive buffer.
    (void)driveSignals;
    (void)lambda;
}


void Foam::NashPanfilovBatched::uploadCUDAInputs() const
{
    cuda_.uploadDriveSignals(driveSignals_.data(), nCells_);
}


void Foam::NashPanfilovBatched::launchCUDAKernel() const
{
    launchNashPanfilovBatchKernel
    (
        cuda_.d_driveSignals,
        cuda_.d_constants,
        nCells_,
        cuda_.d_states,
        cuda_.d_rates,
        cuda_.d_algebraics
    );
}


void Foam::NashPanfilovBatched::downloadCUDATension(scalarField& Ta) const
{
    cuda_.downloadStateRow(Ta.data(), ::Ta, nCells_);
}
#endif


// ************************************************************************* //
