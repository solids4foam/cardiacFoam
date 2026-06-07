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
#include "../LandNiederer/LandNiederer_2017Names.H"
#include "../LandNiederer/LandNiederer_2017.H"

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
    lambdaRate_(num, 0.0)
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
    algebraicValues[AV_lambda_rate] = lambdaRate_[cellI];

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
