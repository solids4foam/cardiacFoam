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

#include "LandNiedererBatched.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "restartStateIO.H"
#include "../LandNiederer/LandNiederer_2017Names.H"
#include "../LandNiederer/LandNiederer_2017.H"

#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(LandNiedererBatched, 0);
    addToRunTimeSelectionTable(activeTensionModel, LandNiedererBatched, dictionary);
}

// * * * * * * * * * * * * * * * * io* hooks  * * * * * * * * * * * * * * * //

const char* const* Foam::LandNiedererBatched::ioStateNames() const
{
    return LandNiedererSTATES_NAMES;
}

const char* const* Foam::LandNiedererBatched::ioConstantNames() const
{
    return LandNiedererCONSTANTS_NAMES;
}

const char* const* Foam::LandNiedererBatched::ioAlgebraicNames() const
{
    return LandNiedererALGEBRAIC_NAMES;
}


bool Foam::LandNiedererBatched::readRestartState(const fvMesh& mesh)
{
    const fileName statePath = restartStateIO::path(mesh, type() + "State");
    if (!isFile(statePath))
    {
        return false;
    }

    syncAllToIO();
    std::ifstream is(statePath.c_str(), std::ios::binary);
    if (!is)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath
            << exit(FatalError);
    }

    restartStateIO::validateHeader
    (
        type(), nStates_ + 1, nCells_, is, statePath
    );

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            ioStates_[cellI][stateI] =
                restartStateIO::readScalar(is, statePath);
            restartStateIO::checkValue(ioStates_[cellI][stateI], statePath);
        }

        prevLambda_[cellI] = restartStateIO::readScalar(is, statePath);
        restartStateIO::checkValue(prevLambda_[cellI], statePath);
    }

    syncStatesFromIO();
    core_.clearTransientSolveData(persistAlgebraics_);
    lambdaRate_ = 0.0;
    return true;
}


void Foam::LandNiedererBatched::writeRestartState(const fvMesh& mesh) const
{
    const fileName statePath = restartStateIO::path(mesh, type() + "State");
    syncAllToIO();
    std::ofstream os(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }

    restartStateIO::writeHeader(type(), nStates_ + 1, nCells_, os);
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            restartStateIO::checkValue(ioStates_[cellI][stateI], statePath);
            restartStateIO::writeScalar(os, ioStates_[cellI][stateI], statePath);
        }
        restartStateIO::checkValue(prevLambda_[cellI], statePath);
        restartStateIO::writeScalar(os, prevLambda_[cellI], statePath);
    }
}


void Foam::LandNiedererBatched::refreshRestartState(const fvMesh& mesh)
{
    const volVectorField& D = mesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = mesh.lookupObject<volVectorField>("f0");
    const volTensorField gradD(fvc::grad(D));
    CellScratch scratch(nStates_, nAlgebraics_, useRushLarsen_);
    BatchedTensionBackend backend(*this);
    const ElectromechanicalSignalProvider& p = provider();

    lambdaRate_ = 0.0;
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        const tensor F(I + gradD[cellI].T());
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            mesh.time().value(),
            p.signal(cellI, driveSignal()),
            mag(F & f0[cellI]),
            scratch
        );
        backend.syncEvaluatedOutputs(cellI, scratch);
        restartTa_[cellI] = backend.activeTensionFromScratch(cellI, scratch);
    }

    ioSynchronized_ = false;
}


bool Foam::LandNiedererBatched::restartTension(scalarField& Ta) const
{
    if (Ta.size() != nCells_)
    {
        return false;
    }

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        Ta[cellI] = restartTa_[cellI];
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LandNiedererBatched::LandNiedererBatched
(
    const dictionary& dict,
    const label num
)
:
    batchedActiveTensionModel(dict, num, NUM_STATES, NUM_ALGEBRAIC),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    prevLambda_(num, 1.0),
    lambdaRate_(num, 0.0),
    restartTa_(num, 0.0)
{
    const word requestedSignal = dict_.lookupOrDefault<word>("couplingSignal", "Cai");

    if (!(requestedSignal == "Cai" || requestedSignal == "cai"))
    {
        FatalErrorInFunction
            << "Unknown LandNiedererBatched 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid option is: Cai."
            << abort(FatalError);
    }

    Info<< nl << "Initialize LandNiedererBatched constants:" << nl;

    scalarField protoStates(NUM_STATES, 0.0);
    scalarField protoRates(NUM_STATES, 0.0);

    // Reuse the 2017 model init hook
    LandNiederer2017initConsts
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
            const word name(ioConstantNames()[k]);
            if (cDict.found(name)) CONSTANTS_[k] = cDict.get<scalar>(name);
        }

        // Re-derive dependent constants
        CONSTANTS_[AC_fPKA_TnI] =
            1.45 - 0.45 * (1.0 - CONSTANTS_[AC_fTnI_PKA]) / (1.0 - CONSTANTS_[AC_fracTnIpo]);
        CONSTANTS_[AC_XSSS] = CONSTANTS_[AC_dr] * 0.5;
        CONSTANTS_[AC_XWSS] = (1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac] * 0.5;
        CONSTANTS_[AC_A] = CONSTANTS_[AC_TOT_A] * CONSTANTS_[AC_dr] /
                           ((1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac] + CONSTANTS_[AC_dr]);
        CONSTANTS_[AC_PKAForceMultiplier] = 1.0 + 0.26 * CONSTANTS_[AC_fMyBPC_PKA];
        CONSTANTS_[AC_k_uw] = 0.026 * CONSTANTS_[AC_nu];
        CONSTANTS_[AC_k_ws] = 0.004 * (1.0 + CONSTANTS_[AC_fMyBPC_PKA] / 2.0) * CONSTANTS_[AC_mu];
        CONSTANTS_[AC_k_wu] = CONSTANTS_[AC_k_uw] * (1.0 / CONSTANTS_[AC_wfrac] - 1.0) - CONSTANTS_[AC_k_ws];
        CONSTANTS_[AC_k_su] = CONSTANTS_[AC_k_ws] * (1.0 / CONSTANTS_[AC_dr] - 1.0) * CONSTANTS_[AC_wfrac];
        CONSTANTS_[AC_cds] = CONSTANTS_[AC_phi] * CONSTANTS_[AC_k_ws] * CONSTANTS_[AC_wfrac] * (1.0 - CONSTANTS_[AC_dr]) / CONSTANTS_[AC_dr];
        CONSTANTS_[AC_cdw] = CONSTANTS_[AC_phi] * CONSTANTS_[AC_k_uw] * (1.0 - CONSTANTS_[AC_wfrac]) / CONSTANTS_[AC_wfrac];
        CONSTANTS_[AC_ktm_block] = CONSTANTS_[AC_ktm_unblock] * std::pow(CONSTANTS_[AC_perm50], CONSTANTS_[AC_nperm]) * 0.5 / (0.5 - CONSTANTS_[AC_XSSS] - CONSTANTS_[AC_XWSS]);
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

void Foam::LandNiedererBatched::preconditionToRestingState
(
    const scalar restingCai
)
{
    const scalar preconditioningTime =
        dict_.lookupOrDefault<scalar>("preconditioningTime", 1000.0); // ms

    if (restingCai < 0)
    {
        FatalErrorInFunction
            << "LandNiedererBatched: resting Ca_i must be non-negative; got "
            << restingCai << " mM. A negative resting calcium is unphysical "
            << "and points to a misconfigured electromechanical signal "
            << "provider." << exit(FatalError);
    }

    // A resting Ca_i at (or below) zero is already the equilibrium of the
    // shipped initial conditions, so preconditioning would be a no-op.
    if (preconditioningTime <= SMALL || restingCai <= SMALL)
    {
        return;
    }

    // Match the scalar LandNiederer model: integrate to resting steady state
    // over a fixed number of substeps rather than a dictionary-configurable
    // step size.
    const label nSteps = 100;
    const scalar step = preconditioningTime/scalar(nSteps);

    CellScratch scratch(nStates_, nAlgebraics_, useRushLarsen_);
    BatchedTensionBackend backend(*this);
    backend.gatherCellState(0, scratch.stateValues);

    for (label stepI = 0; stepI < nSteps; ++stepI)
    {
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            0,
            stepI*step,
            restingCai,
            1.0,
            scratch
        );

        if (useRushLarsen_)
        {
            backend.buildPredictorState(0, step, scratch);
            backend.applyCorrectorState(scratch);
        }
        else
        {
            backend.applyExplicitEulerStep(step, scratch);
        }
    }

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        core_.scatterCellState(cellI, scratch.stateValues);
    }
    core_.clearTransientSolveData(persistAlgebraics_);
    ioSynchronized_ = false;
    prevLambda_ = 1.0;
    lambdaRate_ = 0.0;

    Info<< "    LandNiedererBatched: pre-conditioned " << nCells_
        << " points to resting steady state" << nl
        << "      restingCai = " << restingCai << " mM ("
        << preconditioningTime << " ms integration, "
        << nSteps << " steps)" << nl
        << "      Ca_TRPN   = " << scratch.stateValues[Ca_TRPN] << nl
        << "      TmBlocked = " << scratch.stateValues[TmBlocked] << nl
        << "      XW        = " << scratch.stateValues[XW] << nl
        << "      XS        = " << scratch.stateValues[XS] << nl
        << endl;
}

void Foam::LandNiedererBatched::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField& lambda,
    scalarField& Ta
)
{
    // Compute lambdaRate Field
    for (label i = 0; i < nCells_; ++i)
    {
        if (dt > SMALL)
        {
            lambdaRate_[i] = (lambda[i] - prevLambda_[i]) / dt;
        }
        else
        {
            lambdaRate_[i] = 0.0;
        }
        prevLambda_[i] = lambda[i];
    }

    // Call base class
    batchedActiveTensionModel::calculateTension(t, dt, lambda, Ta);
}


void Foam::LandNiedererBatched::evaluateHotPathStateForCell
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
    // Inject inputs
    algebraicValues[AV_Cai] = driveSignal;
    algebraicValues[AV_lambda] = lambda;
    algebraicValues[AV_lambda_rate] = lambdaRate_[cellI] * 1e-3; // convert s^-1 to ms^-1

    // Compute variables by casting to pointers
    LandNiederer2017computeVariables
    (
        0.0,
        CONSTANTS_.data(),
        rateValues.begin(),
        const_cast<scalar*>(stateValues.begin()),
        algebraicValues.begin()
    );
}

// ************************************************************************* //
