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
#include "addToRunTimeSelectionTable.H"
#include "error.H"

#include "LandNiederer_2017.H"   // self-contained ODE equations

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

defineTypeNameAndDebug(LandNiederer, 0);

addToRunTimeSelectionTable
(
    activeTensionModel,
    LandNiederer,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

LandNiederer::LandNiederer
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

    LandNiederer2017initConsts
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
            const word name(LandNiedererCONSTANTS_NAMES[k]);
            if (cDict.found(name))
            {
                CONSTANTS_[k] = cDict.get<scalar>(name);
            }
        }

        // Re-derive dependent constants after any user override
        // (mirrors the dependency logic in LandNiederer2017initConsts)
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
            const word name(LandNiedererSTATES_NAMES[k]);
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

activeTensionModel::Requirements LandNiederer::requirements() const
{
    Requirements req;
    req.needCai = true;
    req.needVm  = false;
    return req;
}


label LandNiederer::nEqns() const
{
    return NUM_STATES;
}


void LandNiederer::derivatives
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
    ALGEBRAIC_TMP[AV_lambda_rate] = currentLambdaRate_;

    // Compute algebraics and rates (writes into dydt via RATES argument)
    LandNiederer2017computeVariables
    (
        0.0,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<double*>(y.cdata()),  // y is read-only within computeVariables
        ALGEBRAIC_TMP.data()
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
    notImplemented("LandNiederer::jacobian — use an explicit ODE solver (e.g. RK45)");
}


void LandNiederer::solveAtPoint
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
    ALGEBRAIC_i[AV_lambda_rate] = lambda_rate;

    // ------------------------------------------------------------------
    // 4. Advance the ODE system over [currentT_, currentT_ + currentDt_]
    // ------------------------------------------------------------------

    scalarField& STATES_i = STATES_[i];

    scalar tStart = currentT_;
    scalar tEnd   = currentT_ + currentDt_;
    scalar step   = currentDt_;

    odeSolver_->solve(tStart, tEnd, STATES_i, step);

    // ------------------------------------------------------------------
    // 5. Recompute final algebraic values (including Ta)
    // ------------------------------------------------------------------

    LandNiederer2017computeVariables
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
