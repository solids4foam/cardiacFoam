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

#include "../LandNiederer/LandNiederer_2017.H"

#include <cmath>
#include <fstream>

namespace Foam
{

defineTypeNameAndDebug(LandNiedererBatched, 0);
addToRunTimeSelectionTable(activeTensionModel, LandNiedererBatched, dictionary);


const char* const* LandNiedererBatched::ioStateNames() const
{
    return LandNiedererSTATES_NAMES;
}


const char* const* LandNiedererBatched::ioConstantNames() const
{
    return LandNiedererCONSTANTS_NAMES;
}


const char* const* LandNiedererBatched::ioAlgebraicNames() const
{
    return LandNiedererALGEBRAIC_NAMES;
}


void LandNiedererBatched::updateDerivedConstants()
{
    CONSTANTS_[AC_k_ws] = 0.004 * CONSTANTS_[AC_mu];
    CONSTANTS_[AC_k_uw] = 0.026 * CONSTANTS_[AC_nu];

    CONSTANTS_[AC_cdw] =
        CONSTANTS_[AC_phi]
      * CONSTANTS_[AC_k_uw]
      * (1.0 - CONSTANTS_[AC_dr])
      * (1.0 - CONSTANTS_[AC_wfrac])
      / ((1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac]);

    CONSTANTS_[AC_cds] =
        CONSTANTS_[AC_phi]
      * CONSTANTS_[AC_k_ws]
      * (1.0 - CONSTANTS_[AC_dr])
      * CONSTANTS_[AC_wfrac]
      / CONSTANTS_[AC_dr];

    CONSTANTS_[AC_k_wu] =
        CONSTANTS_[AC_k_uw]
      * (1.0 / CONSTANTS_[AC_wfrac] - 1.0)
      - CONSTANTS_[AC_k_ws];

    CONSTANTS_[AC_k_su] =
        CONSTANTS_[AC_k_ws]
      * (1.0 / CONSTANTS_[AC_dr] - 1.0)
      * CONSTANTS_[AC_wfrac];

    CONSTANTS_[AC_A] =
        (0.25 * CONSTANTS_[AC_TOT_A])
      / ((1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac]
         + CONSTANTS_[AC_dr])
      * (CONSTANTS_[AC_dr] / 0.25);

    CONSTANTS_[AC_XSSS] = CONSTANTS_[AC_dr] * 0.5;
    CONSTANTS_[AC_XWSS] =
        (1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac] * 0.5;

    CONSTANTS_[AC_ktm_block] =
        CONSTANTS_[AC_ktm_unblock]
      * std::pow(CONSTANTS_[AC_perm50], CONSTANTS_[AC_nperm])
      * 0.5
      / (0.5 - CONSTANTS_[AC_XSSS] - CONSTANTS_[AC_XWSS]);
}


LandNiedererBatched::LandNiedererBatched
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
            << requestedSignal << nl << "Valid option is: Cai."
            << exit(FatalError);
    }

    scalarField protoStates(NUM_STATES, 0.0);
    scalarField protoRates(NUM_STATES, 0.0);
    LandNiederer2017initConsts
    (
        CONSTANTS_.data(),
        protoRates.data(),
        protoStates.data()
    );

    if (dict_.found("constants"))
    {
        const dictionary& constants = dict_.subDict("constants");
        for (label constantI = 0; constantI < NUM_CONSTANTS; ++constantI)
        {
            const word name(LandNiedererCONSTANTS_NAMES[constantI]);
            if (constants.found(name))
            {
                CONSTANTS_[constantI] = constants.get<scalar>(name);
            }
        }
        updateDerivedConstants();
    }

    if (dict_.found("initialStates"))
    {
        const dictionary& initialStates = dict_.subDict("initialStates");
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            const word name(LandNiedererSTATES_NAMES[stateI]);
            if (initialStates.found(name))
            {
                protoStates[stateI] = initialStates.get<scalar>(name);
            }
        }
    }

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = protoStates[stateI];
        }
    }
}


void LandNiedererBatched::preconditionToRestingState(const scalarField& restingCai)
{
    const scalar preconditioningTime =
        dict_.lookupOrDefault<scalar>("preconditioningTime", 1000.0);
    if (restingCai.size() != nCells_)
    {
        FatalErrorInFunction
            << "LandNiedererBatched received " << restingCai.size()
            << " resting Ca_i values for " << nCells_ << " cells."
            << exit(FatalError);
    }
    if (min(restingCai) < 0.0)
    {
        FatalErrorInFunction
            << "LandNiedererBatched requires a non-negative resting Ca_i, got "
            << min(restingCai) << " mM." << exit(FatalError);
    }
    if (preconditioningTime <= SMALL)
    {
        return;
    }

    const label nSteps = 100;
    const scalar dt = preconditioningTime / scalar(nSteps);
    const bool uniform = (max(restingCai) - min(restingCai)) <= SMALL;
    CellScratch scratch(nStates_, nAlgebraics_);
    BatchedTensionBackend backend(*this);

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        if (!(uniform && cellI > 0))
        {
            backend.gatherCellState(cellI, scratch.stateValues);

            for (label stepI = 0; stepI < nSteps; ++stepI)
            {
                scratch.resetPrimary();
                backend.evaluateScratchAtTime
                (
                    cellI,
                    stepI * dt,
                    scaledDriveSignal(restingCai[cellI]),
                    1.0,
                    scratch
                );

                backend.advanceSubstep(cellI, dt, scratch);
            }
        }

        core_.scatterCellState(cellI, scratch.stateValues);
    }

    core_.clearTransientSolveData(persistAlgebraics_);
    ioSynchronized_ = false;
    prevLambda_ = 1.0;
    lambdaRate_ = 0.0;
}


void LandNiedererBatched::calculateTension
(
    const scalar t,
    const scalar dt,
    const scalarField& lambda,
    scalarField& Ta
)
{
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        lambdaRate_[cellI] =
            dt > SMALL ? (lambda[cellI] - prevLambda_[cellI])/dt : 0.0;
        prevLambda_[cellI] = lambda[cellI];
    }

    batchedActiveTensionModel::calculateTension(t, dt, lambda, Ta);
}


void LandNiedererBatched::evaluateHotPathStateForCell
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
    algebraicValues[AV_Cai] = driveSignal;
    algebraicValues[AV_lambda] = lambda;
    algebraicValues[AV_lambda_rate] = lambdaRate_[cellI] * 1.0e-3;

    LandNiederer2017computeVariables
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.begin(),
        const_cast<scalar*>(stateValues.begin()),
        algebraicValues.begin()
    );
}


bool LandNiedererBatched::readRestartState(const fvMesh& mesh)
{
    const fileName statePath = restartStateIO::path(mesh, type() + "State");
    if (!isFile(statePath))
    {
        return false;
    }

    syncAllToIO();
    std::ifstream input(statePath.c_str(), std::ios::binary);
    if (!input)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath << exit(FatalError);
    }
    restartStateIO::validateHeader
    (
        type(), nStates_ + 1, nCells_, input, statePath
    );

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            ioStates_[cellI][stateI] = restartStateIO::readScalar(input, statePath);
            restartStateIO::checkValue(ioStates_[cellI][stateI], statePath);
        }
        prevLambda_[cellI] = restartStateIO::readScalar(input, statePath);
        restartStateIO::checkValue(prevLambda_[cellI], statePath);
    }

    syncStatesFromIO();
    core_.clearTransientSolveData(persistAlgebraics_);
    lambdaRate_ = 0.0;
    return true;
}


void LandNiedererBatched::writeRestartState(const fvMesh& mesh) const
{
    const fileName statePath = restartStateIO::path(mesh, type() + "State");
    syncAllToIO();
    std::ofstream output(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath << exit(FatalError);
    }
    restartStateIO::writeHeader(type(), nStates_ + 1, nCells_, output);

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < nStates_; ++stateI)
        {
            restartStateIO::checkValue(ioStates_[cellI][stateI], statePath);
            restartStateIO::writeScalar(output, ioStates_[cellI][stateI], statePath);
        }
        restartStateIO::checkValue(prevLambda_[cellI], statePath);
        restartStateIO::writeScalar(output, prevLambda_[cellI], statePath);
    }
}


void LandNiedererBatched::refreshRestartState(const fvMesh& mesh)
{
    const volVectorField& D = mesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = mesh.lookupObject<volVectorField>("f0");
    const volTensorField gradD(fvc::grad(D));
    CellScratch scratch(nStates_, nAlgebraics_);
    BatchedTensionBackend backend(*this);
    lambdaRate_ = 0.0;
    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        const tensor F(I + gradD[cellI].T());
        backend.gatherCellState(cellI, scratch.stateValues);
        scratch.resetPrimary();
        backend.evaluateScratchAtTime
        (
            cellI,
            mesh.time().value()*timeScaleFactor(),
            coupledDriveSignal(cellI),
            mag(F & f0[cellI]),
            scratch
        );
        backend.syncEvaluatedOutputs(cellI, scratch);
        restartTa_[cellI] = backend.activeTensionFromScratch(cellI, scratch);
    }

    ioSynchronized_ = false;
}


bool LandNiedererBatched::restartTension(scalarField& Ta) const
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

} // End namespace Foam

// ************************************************************************* //
