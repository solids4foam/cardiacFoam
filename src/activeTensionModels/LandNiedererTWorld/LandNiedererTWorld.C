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

#include "LandNiedererTWorld.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "fvcGrad.H"
#include "restartStateIO.H"

#include "LandNiedererTWorld_2025.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

defineTypeNameAndDebug(LandNiedererTWorld, 0);

addToRunTimeSelectionTable
(
    activeTensionModel,
    LandNiedererTWorld,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const char* const* LandNiedererTWorld::ioConstantNames() const
{
    return LandNiedererTWorldCONSTANTS_NAMES;
}

const char* const* LandNiedererTWorld::ioStateNames() const
{
    return LandNiedererTWorldSTATES_NAMES;
}

const char* const* LandNiedererTWorld::ioAlgebraicNames() const
{
    return LandNiedererTWorldALGEBRAIC_NAMES;
}


bool LandNiedererTWorld::readRestartState(const fvMesh& mesh)
{
    const fileName statePath =
        restartStateIO::path(mesh, "LandNiedererTWorldState");

    if (!isFile(statePath))
    {
        return false;
    }

    std::ifstream is(statePath.c_str(), std::ios::binary);
    if (!is)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath
            << exit(FatalError);
    }
    restartStateIO::validateHeader
    (
        "LandNiedererTWorld", NUM_STATES + 1, STATES_.size(), is, statePath
    );

    forAll(STATES_, integrationPtI)
    {
        forAll(STATES_[integrationPtI], stateI)
        {
            STATES_[integrationPtI][stateI] =
                restartStateIO::readScalar(is, statePath);
            restartStateIO::checkValue
            (
                STATES_[integrationPtI][stateI], statePath
            );
        }

        prevLambda_[integrationPtI] =
            restartStateIO::readScalar(is, statePath);
        restartStateIO::checkValue(prevLambda_[integrationPtI], statePath);
    }

    const volVectorField& D = mesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = mesh.lookupObject<volVectorField>("f0");
    const volTensorField gradD(fvc::grad(D));
    const scalar t = mesh.time().value()*1000.0;

    forAll(STATES_, integrationPtI)
    {
        const tensor F(I + gradD[integrationPtI].T());
        const scalar lambda = mag(F & f0[integrationPtI]);
        const scalar driveVal = provider().signal
        (
            integrationPtI, CouplingSignal::CAI
        );

        currentDriveSignal_ = driveVal;
        currentLambda_ = lambda;
        currentLambdaRate_ = 0.0;
        ALGEBRAIC_[integrationPtI][AV_Cai] = driveVal;
        ALGEBRAIC_[integrationPtI][AV_lambda] = lambda;
        ALGEBRAIC_[integrationPtI][AV_lambda_rate] = 0.0;

        LandNiedererTWorld2025computeVariables
        (
            t,
            CONSTANTS_.data(),
            RATES_[integrationPtI].data(),
            STATES_[integrationPtI].data(),
            ALGEBRAIC_[integrationPtI].data()
        );
    }

    return true;
}


void LandNiedererTWorld::writeRestartState(const fvMesh& mesh) const
{
    const fileName statePath =
        restartStateIO::path(mesh, "LandNiedererTWorldState");
    std::ofstream os(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }
    restartStateIO::writeHeader
    (
        "LandNiedererTWorld", NUM_STATES + 1, STATES_.size(), os
    );

    forAll(STATES_, integrationPtI)
    {
        forAll(STATES_[integrationPtI], stateI)
        {
            restartStateIO::checkValue
            (
                STATES_[integrationPtI][stateI], statePath
            );
            restartStateIO::writeScalar
            (
                os, STATES_[integrationPtI][stateI], statePath
            );
        }
        restartStateIO::checkValue(prevLambda_[integrationPtI], statePath);
        restartStateIO::writeScalar(os, prevLambda_[integrationPtI], statePath);
    }
}


bool LandNiedererTWorld::restartTension(scalarField& Ta) const
{
    if (Ta.size() != ALGEBRAIC_.size())
    {
        FatalErrorInFunction
            << "Restart tension size mismatch"
            << exit(FatalError);
    }

    forAll(Ta, integrationPtI)
    {
        Ta[integrationPtI] = ALGEBRAIC_[integrationPtI][AV_Ta];
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

LandNiedererTWorld::LandNiedererTWorld
(
    const dictionary& dict,
    const label nIntegrationPoints
)
:
    activeTensionModel(dict, nIntegrationPoints),
    STATES_(nIntegrationPoints),
    ALGEBRAIC_(nIntegrationPoints),
    RATES_(nIntegrationPoints),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    prevLambda_(nIntegrationPoints, 1.0),
    currentLambda_(1.0),
    currentLambdaRate_(0.0)
{
    // Initialise shared constants and prototype state values
    scalarField protoStates(NUM_STATES, 0.0);
    scalarField protoRates(NUM_STATES, 0.0);

    LandNiedererTWorld2025initConsts
    (
        CONSTANTS_.data(),
        protoRates.data(),
        protoStates.data()
    );

    // Allow the user to override individual constants from the dictionary
    if (dict.found("constants"))
    {
        const dictionary& cDict = dict.subDict("constants");

        for (label k = 0; k < NUM_CONSTANTS; ++k)
        {
            const word name(LandNiedererTWorldCONSTANTS_NAMES[k]);
            if (cDict.found(name))
            {
                CONSTANTS_[k] = cDict.get<scalar>(name);
            }
        }

        // Re-derive dependent constants after any user override
        // (mirrors the dependency logic in LandNiedererTWorld2025initConsts)
        CONSTANTS_[AC_fPKA_TnI] =
            1.45
          - 0.45 * (1.0 - CONSTANTS_[AC_fTnI_PKA])
                 / (1.0 - CONSTANTS_[AC_fracTnIpo]);

        CONSTANTS_[AC_XSSS] = CONSTANTS_[AC_dr] * 0.5;
        CONSTANTS_[AC_XWSS] =
            (1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac] * 0.5;

        CONSTANTS_[AC_A] =
            CONSTANTS_[AC_TOT_A] * CONSTANTS_[AC_dr]
          / (  (1.0 - CONSTANTS_[AC_dr]) * CONSTANTS_[AC_wfrac]
             + CONSTANTS_[AC_dr]);

        CONSTANTS_[AC_PKAForceMultiplier] =
            1.0 + 0.26 * CONSTANTS_[AC_fMyBPC_PKA];

        CONSTANTS_[AC_k_uw] = 0.026 * CONSTANTS_[AC_nu];

        CONSTANTS_[AC_k_ws] =
            0.004
          * (1.0 + CONSTANTS_[AC_fMyBPC_PKA] / 2.0)
          * CONSTANTS_[AC_mu];

        CONSTANTS_[AC_k_wu] =
            CONSTANTS_[AC_k_uw] * (1.0 / CONSTANTS_[AC_wfrac] - 1.0)
          - CONSTANTS_[AC_k_ws];

        CONSTANTS_[AC_k_su] =
            CONSTANTS_[AC_k_ws]
          * (1.0 / CONSTANTS_[AC_dr] - 1.0)
          * CONSTANTS_[AC_wfrac];

        CONSTANTS_[AC_cds] =
            CONSTANTS_[AC_phi] * CONSTANTS_[AC_k_ws]
          * CONSTANTS_[AC_wfrac] * (1.0 - CONSTANTS_[AC_dr])
          / CONSTANTS_[AC_dr];

        CONSTANTS_[AC_cdw] =
            CONSTANTS_[AC_phi] * CONSTANTS_[AC_k_uw]
          * (1.0 - CONSTANTS_[AC_wfrac]) / CONSTANTS_[AC_wfrac];

        CONSTANTS_[AC_ktm_block] =
            CONSTANTS_[AC_ktm_unblock]
          * std::pow(CONSTANTS_[AC_perm50], CONSTANTS_[AC_nperm])
          * 0.5
          / (0.5 - CONSTANTS_[AC_XSSS] - CONSTANTS_[AC_XWSS]);
    }

    // Allow the user to override initial state values from the dictionary
    if (dict.found("initialStates"))
    {
        const dictionary& sDict = dict.subDict("initialStates");

        for (label k = 0; k < NUM_STATES; ++k)
        {
            const word name(LandNiedererTWorldSTATES_NAMES[k]);
            if (sDict.found(name))
            {
                protoStates[k] = sDict.get<scalar>(name);
            }
        }
    }

    // Allocate and initialise per-integration-point storage
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(protoStates));
        ALGEBRAIC_.set(i, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(i, new scalarField(NUM_STATES, 0.0));
    }

    // Build ODE solver
    const dictionary& solverDict =
        dict.found("ODESolver")
      ? dict.subDict("ODESolver")
      : dict;                                  // fallback to outer dict

    odeSolver_ = ODESolver::New(*this, solverDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void LandNiedererTWorld::preconditionToRestingState(const scalar restingCai)
{
    // Run the ODE for preconditioningTime ms at constant Ca_i=restingCai,
    // lambda=1, lambda_rate=0 to drive all states to resting equilibrium.
    //
    // The TWorld contraction initial conditions (Ca_TRPN=0, XS=0) are the
    // equilibrium at Ca_i=0.  At the ionic model's actual resting Ca_i
    // (e.g. TNNP: 0.2 µM, TWorld: 0.097 µM, ToRORd: 0.075 µM) the true
    // resting XS is non-zero, so starting from 0 produces a spurious global
    // Ta transient in all cells simultaneously.
    const scalar preconditioningTime =
        dict_.lookupOrDefault<scalar>("preconditioningTime", 1000.0); // ms

    if (restingCai < 0)
    {
        FatalErrorInFunction
            << "LandNiedererTWorld: resting Ca_i must be non-negative; got "
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

    currentDriveSignal_ = restingCai;
    currentLambda_      = 1.0;
    currentLambdaRate_  = 0.0;

    // Run a single-point ODE integration; all cells share the same resting IC.
    scalarField precondStates(STATES_[0]);
    scalar step = preconditioningTime / 100.0;
    odeSolver_->solve(0.0, preconditioningTime, precondStates, step);

    forAll(STATES_, i)
    {
        STATES_[i] = precondStates;
    }

    Info<< "    LandNiedererTWorld: pre-conditioned " << STATES_.size()
        << " points to resting steady state" << nl
        << "      restingCai = " << restingCai << " mM ("
        << preconditioningTime << " ms integration)" << nl
        << "      Ca_TRPN   = " << precondStates[Ca_TRPN] << nl
        << "      TmBlocked = " << precondStates[TmBlocked] << nl
        << "      XW        = " << precondStates[XW] << nl
        << "      XS        = " << precondStates[XS] << nl
        << endl;
}


activeTensionModel::Requirements LandNiedererTWorld::requirements() const
{
    Requirements req;
    req.needCai     = true;
    req.needVm      = false;
    req.needsLambda = true;
    return req;
}


label LandNiedererTWorld::nEqns() const
{
    return NUM_STATES;
}


void LandNiedererTWorld::derivatives
(
    const scalar /*t*/,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Temporary algebraic storage for this sub-step
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    // Inject drive signal and stretch into temporary algebraics
    ALGEBRAIC_TMP[AV_Cai]         = currentDriveSignal_;  // [Ca2+]_i, mM
    ALGEBRAIC_TMP[AV_lambda]      = currentLambda_;
    ALGEBRAIC_TMP[AV_lambda_rate] = currentLambdaRate_ * 1e-3; // convert s^-1 to ms^-1

    // Compute algebraics and rates (writes into dydt via RATES argument)
    LandNiedererTWorld2025computeVariables
    (
        0.0,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<double*>(y.cdata()),  // y is read-only within computeVariables
        ALGEBRAIC_TMP.data()
    );
}

void LandNiedererTWorld::jacobian
(
    const scalar,
    const scalarField&,
    scalarField&,
    scalarSquareMatrix&
) const
{
    notImplemented
    (
        "LandNiedererTWorld::jacobian — use an explicit ODE solver (e.g. RK45)"
    );
}


void LandNiedererTWorld::solveAtPoint
(
    label i,
    scalar driveVal,
    scalar lambda,
    scalar& Ta
)
{
    // ------------------------------------------------------------------
    // 1. Compute lambda_rate from previous stored value
    // ------------------------------------------------------------------

    scalar lambda_rate = 0.0;

    if (currentDt_ > SMALL)
    {
        lambda_rate = (lambda - prevLambda_[i]) / currentDt_;

        // Cap to a physiologically plausible range (Land 2017 calibrated at
        // ~1 ms steps; rates beyond ±20 s^-1 indicate a predictor overshoot
        // or first-activation transient, not real sarcomere dynamics).
        lambda_rate = max(min(lambda_rate, scalar(20.0)), scalar(-20.0));
    }

    prevLambda_[i] = lambda;

    // ------------------------------------------------------------------
    // 2. Store per-point quantities for use in derivatives() sub-steps
    // ------------------------------------------------------------------

    currentDriveSignal_ = driveVal;   // base class member: [Ca2+]_i in mM
    currentLambda_      = lambda;
    currentLambdaRate_  = lambda_rate;

    // ------------------------------------------------------------------
    // 3. Inject inputs into algebraic storage
    // ------------------------------------------------------------------

    scalarField& ALGEBRAIC_i = ALGEBRAIC_[i];
    ALGEBRAIC_i[AV_Cai]         = driveVal;
    ALGEBRAIC_i[AV_lambda]      = lambda;
    ALGEBRAIC_i[AV_lambda_rate] = lambda_rate * 1e-3; // convert s^-1 to ms^-1

    // ------------------------------------------------------------------
    // 4. Advance the ODE system over [currentT_, currentT_ + currentDt_]
    // ------------------------------------------------------------------

    scalarField& STATES_i = STATES_[i];

    // The TWorld contraction subsystem operates natively in ms, so we must
    // scale the integration time limits by 1000 to solve the ODE in ms.
    const scalar tStart = currentT_ * 1000.0;
    const scalar tEnd   = (currentT_ + currentDt_) * 1000.0;
    scalar step         = currentDt_ * 1000.0;

    odeSolver_->solve(tStart, tEnd, STATES_i, step);

    // ------------------------------------------------------------------
    // 5. Recompute final algebraic values (including Ta)
    // ------------------------------------------------------------------

    LandNiedererTWorld2025computeVariables
    (
        tEnd,
        CONSTANTS_.data(),
        RATES_[i].data(),
        STATES_i.data(),
        ALGEBRAIC_i.data()
    );

    // ------------------------------------------------------------------
    // 6. Extract active tension (algebraic output)
    // ------------------------------------------------------------------

    Ta = ALGEBRAIC_i[AV_Ta];
}

} // End namespace Foam

// ************************************************************************* //
