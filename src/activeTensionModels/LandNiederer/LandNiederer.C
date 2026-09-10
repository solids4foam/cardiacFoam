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

#include "LandNiederer.H"
#include "LandNiederer_2017.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "fvcGrad.H"
#include "restartStateIO.H"

#include <fstream>

namespace Foam
{

defineTypeNameAndDebug(LandNiederer, 0);
addToRunTimeSelectionTable(activeTensionModel, LandNiederer, dictionary);


const char* const* LandNiederer::ioConstantNames() const
{
    return LandNiedererCONSTANTS_NAMES;
}


const char* const* LandNiederer::ioStateNames() const
{
    return LandNiedererSTATES_NAMES;
}


const char* const* LandNiederer::ioAlgebraicNames() const
{
    return LandNiedererALGEBRAIC_NAMES;
}


void LandNiederer::updateDerivedConstants()
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


LandNiederer::LandNiederer
(
    const dictionary& dict,
    const label nIntegrationPoints
)
:
    activeTensionModel(dict, nIntegrationPoints),
    odeSolver_(),
    STATES_(nIntegrationPoints),
    ALGEBRAIC_(nIntegrationPoints),
    RATES_(nIntegrationPoints),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    prevLambda_(nIntegrationPoints, 1.0),
    currentLambda_(1.0),
    currentLambdaRate_(0.0)
{
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

    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI, new scalarField(protoStates));
        ALGEBRAIC_.set
        (
            integrationPtI,
            new scalarField(NUM_ALGEBRAIC, 0.0)
        );
        RATES_.set(integrationPtI, new scalarField(NUM_STATES, 0.0));
    }

    const dictionary& solverDict =
        dict_.found("ODESolver") ? dict_.subDict("ODESolver") : dict_;
    odeSolver_ = ODESolver::New(*this, solverDict);
}


activeTensionModel::Requirements LandNiederer::requirements() const
{
    Requirements requirements;
    requirements.needCai = true;
    requirements.needsLambda = true;
    return requirements;
}


void LandNiederer::preconditionToRestingState(const scalarField& restingCai)
{
    const scalar preconditioningTime =
        dict_.lookupOrDefault<scalar>("preconditioningTime", 1000.0);

    if (restingCai.size() != STATES_.size())
    {
        FatalErrorInFunction
            << "LandNiederer received " << restingCai.size()
            << " resting Ca_i values for " << STATES_.size()
            << " integration points." << exit(FatalError);
    }

    if (min(restingCai) < 0.0)
    {
        FatalErrorInFunction
            << "LandNiederer requires a non-negative resting Ca_i, got "
            << min(restingCai) << " mM." << exit(FatalError);
    }

    if (preconditioningTime <= SMALL)
    {
        return;
    }

    currentLambda_ = 1.0;
    currentLambdaRate_ = 0.0;

    const bool uniform = (max(restingCai) - min(restingCai)) <= SMALL;
    scalarField preconditionedStates(STATES_[0]);

    forAll(STATES_, integrationPtI)
    {
        if (!(uniform && integrationPtI > 0))
        {
            currentDriveSignal_ = scaledDriveSignal(restingCai[integrationPtI]);
            preconditionedStates = STATES_[integrationPtI];
            scalar step = preconditioningTime / 100.0;
            odeSolver_->solve(0.0, preconditioningTime, preconditionedStates, step);
        }

        STATES_[integrationPtI] = preconditionedStates;
        prevLambda_[integrationPtI] = 1.0;
    }
}


void LandNiederer::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraic(NUM_ALGEBRAIC, 0.0);
    algebraic[AV_Cai] = currentDriveSignal_;
    algebraic[AV_lambda] = currentLambda_;
    algebraic[AV_lambda_rate] = currentLambdaRate_ * 1.0e-3;

    LandNiederer2017computeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        algebraic.data()
    );
}


void LandNiederer::jacobian
(
    const scalar,
    const scalarField&,
    scalarField&,
    scalarSquareMatrix&
) const
{
    notImplemented
    (
        "LandNiederer::jacobian — use an explicit ODE solver (for example RKF45)"
    );
}


void LandNiederer::solveAtPoint
(
    const label i,
    const scalar driveVal,
    const scalar lambda,
    scalar& Ta
)
{
    if (driveVal < 0.0)
    {
        FatalErrorInFunction
            << "LandNiederer requires a non-negative Ca_i, got " << driveVal
            << " uM at integration point " << i << '.' << exit(FatalError);
    }

    scalar lambdaRate = 0.0;
    if (currentDt_ > SMALL)
    {
        lambdaRate = (lambda - prevLambda_[i])/currentDt_;
    }
    prevLambda_[i] = lambda;

    currentDriveSignal_ = driveVal;
    currentLambda_ = lambda;
    currentLambdaRate_ = lambdaRate;

    scalarField& states = STATES_[i];
    scalarField& algebraic = ALGEBRAIC_[i];
    scalarField& rates = RATES_[i];

    algebraic[AV_Cai] = driveVal;
    algebraic[AV_lambda] = lambda;
    algebraic[AV_lambda_rate] = lambdaRate * 1.0e-3;

    const scalar tStart = scaledTime(currentT_);
    const scalar tEnd = scaledTime(currentT_ + currentDt_);
    scalar step = scaledTime(currentDt_);
    odeSolver_->solve(tStart, tEnd, states, step);

    LandNiederer2017computeVariables
    (
        tEnd,
        CONSTANTS_.data(),
        rates.data(),
        states.data(),
        algebraic.data()
    );

    // Expose active tension only.
    Ta = algebraic[AV_Ta];
}


bool LandNiederer::readRestartState(const fvMesh& mesh)
{
    const fileName statePath = restartStateIO::path(mesh, "LandNiedererState");
    if (!isFile(statePath))
    {
        return false;
    }

    std::ifstream input(statePath.c_str(), std::ios::binary);
    if (!input)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath << exit(FatalError);
    }
    restartStateIO::validateHeader
    (
        "LandNiederer", NUM_STATES + 1, STATES_.size(), input, statePath
    );

    forAll(STATES_, integrationPtI)
    {
        forAll(STATES_[integrationPtI], stateI)
        {
            STATES_[integrationPtI][stateI] =
                restartStateIO::readScalar(input, statePath);
            restartStateIO::checkValue(STATES_[integrationPtI][stateI], statePath);
        }
        prevLambda_[integrationPtI] = restartStateIO::readScalar(input, statePath);
        restartStateIO::checkValue(prevLambda_[integrationPtI], statePath);
    }

    const volVectorField& D = mesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = mesh.lookupObject<volVectorField>("f0");
    const volTensorField gradD(fvc::grad(D));
    const scalar t = scaledTime(mesh.time().value());

    forAll(STATES_, integrationPtI)
    {
        const tensor F(I + gradD[integrationPtI].T());
        const scalar lambda = mag(F & f0[integrationPtI]);
        const scalar cai = coupledDriveSignal(integrationPtI);
        scalarField& algebraic = ALGEBRAIC_[integrationPtI];

        algebraic[AV_Cai] = cai;
        algebraic[AV_lambda] = lambda;
        algebraic[AV_lambda_rate] = 0.0;
        LandNiederer2017computeVariables
        (
            t,
            CONSTANTS_.data(),
            RATES_[integrationPtI].data(),
            STATES_[integrationPtI].data(),
            algebraic.data()
        );
    }

    return true;
}


void LandNiederer::writeRestartState(const fvMesh& mesh) const
{
    const fileName statePath = restartStateIO::path(mesh, "LandNiedererState");
    std::ofstream output(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath << exit(FatalError);
    }
    restartStateIO::writeHeader
    (
        "LandNiederer", NUM_STATES + 1, STATES_.size(), output
    );

    forAll(STATES_, integrationPtI)
    {
        forAll(STATES_[integrationPtI], stateI)
        {
            restartStateIO::checkValue(STATES_[integrationPtI][stateI], statePath);
            restartStateIO::writeScalar
            (
                output, STATES_[integrationPtI][stateI], statePath
            );
        }
        restartStateIO::checkValue(prevLambda_[integrationPtI], statePath);
        restartStateIO::writeScalar(output, prevLambda_[integrationPtI], statePath);
    }
}


bool LandNiederer::restartTension(scalarField& Ta) const
{
    if (Ta.size() != ALGEBRAIC_.size())
    {
        return false;
    }

    forAll(Ta, integrationPtI)
    {
        Ta[integrationPtI] = ALGEBRAIC_[integrationPtI][AV_Ta];
    }
    return true;
}

} // End namespace Foam

// ************************************************************************* //
