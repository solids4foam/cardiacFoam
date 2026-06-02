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

Description
    State-major (SoA) batched evaluator for the reduced Bueno-Orovio 2008
    ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "BuenoOrovio_2008Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

CARDIAC_HOST_DEVICE inline double bo_smoothHeaviside
(
    const double x,
    const double width = 1.0e-4
)
{
    return 0.5 * (1.0 + std::tanh(x / width));
}


namespace Foam
{

enum BO_BATCH_SUPPORT_INDEX
{
    BO_BATCH_SUPPORT_tau_v,
    BO_BATCH_SUPPORT_gInf_v,
    BO_BATCH_SUPPORT_tau_w,
    BO_BATCH_SUPPORT_gInf_w,
    BO_BATCH_SUPPORT_tau_s,
    BO_BATCH_SUPPORT_gInf_s,
    BO_BATCH_SUPPORT_Iion,

    NUM_BO_BATCH_SUPPORT
};

static_assert(NUM_BO_BATCH_SUPPORT == 7, "BuenoOrovio compact support layout drift");


CARDIAC_HOST_DEVICE inline void BuenoOrovioComputeVariablesBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVm,
    const StimulusProtocolPOD& stimulus
)
{
    const double Istim = computeStimulusPOD(t, stimulus);

    const double c_uO         = CONSTANTS[uO];
    const double c_uU         = CONSTANTS[uU];
    const double c_thetaV     = CONSTANTS[thetaV];
    const double c_thetaW     = CONSTANTS[thetaW];
    const double c_thetaVm    = CONSTANTS[thetaVMinus];
    const double c_thetaO     = CONSTANTS[thetaO];
    const double c_tauV1m     = CONSTANTS[tauV1Minus];
    const double c_tauV2m     = CONSTANTS[tauV2Minus];
    const double c_tauVPlus   = CONSTANTS[tauVPlus];
    const double c_tauW1m     = CONSTANTS[tauW1Minus];
    const double c_tauW2m     = CONSTANTS[tauW2Minus];
    const double c_kWm        = CONSTANTS[kWMinus];
    const double c_uWm        = CONSTANTS[uWMinus];
    const double c_tauWPlus   = CONSTANTS[tauWPlus];
    const double c_tauFi      = CONSTANTS[tauFi];
    const double c_tauO1      = CONSTANTS[tauO1];
    const double c_tauO2      = CONSTANTS[tauO2];
    const double c_tauSo1     = CONSTANTS[tauSo1];
    const double c_tauSo2     = CONSTANTS[tauSo2];
    const double c_kSo        = CONSTANTS[kSo];
    const double c_uSo        = CONSTANTS[uSo];
    const double c_tauS1      = CONSTANTS[tauS1];
    const double c_tauS2      = CONSTANTS[tauS2];
    const double c_kS         = CONSTANTS[kS];
    const double c_uS         = CONSTANTS[uS];
    const double c_tauSi      = CONSTANTS[tauSi];
    const double c_tauWInfty  = CONSTANTS[tauWInfty];
    const double c_wInftyStar = CONSTANTS[wInftyStar];

    const int U_BASE = u * N;
    const int V_BASE = v * N;
    const int W_BASE = w * N;
    const int S_BASE = s * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double Su = STATES_SoA[U_BASE + cellI];
        const double Sv = STATES_SoA[V_BASE + cellI];
        const double Sw = STATES_SoA[W_BASE + cellI];
        const double Ss = STATES_SoA[S_BASE + cellI];

        const double hV  = bo_smoothHeaviside(Su - c_thetaV);
        const double hW  = bo_smoothHeaviside(Su - c_thetaW);
        const double hVm = bo_smoothHeaviside(Su - c_thetaVm);
        const double hO  = bo_smoothHeaviside(Su - c_thetaO);

        const double tauSoLocal =
            c_tauSo1 + (c_tauSo2 - c_tauSo1)
            * (1.0 + std::tanh(c_kSo * (Su - c_uSo))) * 0.5;

        const double tauOLocal = (1.0 - hO) * c_tauO1 + hO * c_tauO2;

        const double tauVMinusLocal = (1.0 - hVm) * c_tauV1m + hVm * c_tauV2m;

        const double tauWMinusLocal =
            c_tauW1m + (c_tauW2m - c_tauW1m)
            * (1.0 + std::tanh(c_kWm * (Su - c_uWm))) * 0.5;

        const double tauSLocal = (1.0 - hW) * c_tauS1 + hW * c_tauS2;

        const double Jfi =
            -Sv * hV * (Su - c_thetaV) * (c_uU - Su) / c_tauFi;

        const double Jso =
            (Su - c_uO) * (1.0 - hW) / tauOLocal
          + hW / tauSoLocal;

        const double Jsi = -hW * Sw * Ss / c_tauSi;

        const double Jion = Jfi + Jso + Jsi;

        const double vInftyLocal = bo_smoothHeaviside(c_thetaVm - Su);

        const double wInftyLocal =
            (1.0 - hO) * (1.0 - Su / c_tauWInfty)
          + hO * c_wInftyStar;

        const double invTauV =
            (1.0 - hV)/tauVMinusLocal + hV/c_tauVPlus;
        const double tauVEffective = 1.0/invTauV;
        const double vSteady =
            (1.0 - hV)*vInftyLocal/(tauVMinusLocal*invTauV);

        const double invTauW =
            (1.0 - hW)/tauWMinusLocal + hW/c_tauWPlus;
        const double tauWEffective = 1.0/invTauW;
        const double wSteady =
            (1.0 - hW)*wInftyLocal/(tauWMinusLocal*invTauW);

        const double sSteady =
            (1.0 + std::tanh(c_kS * (Su - c_uS))) * 0.5;

        RATES_SoA[V_BASE + cellI] = (vSteady - Sv)/tauVEffective;
        RATES_SoA[W_BASE + cellI] = (wSteady - Sw)/tauWEffective;
        RATES_SoA[S_BASE + cellI] = (sSteady - Ss)/tauSLocal;

        if (solveVm)
        {
            RATES_SoA[U_BASE + cellI] = -Jion - Istim;
        }

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[BO_BATCH_SUPPORT_tau_v  * N + cellI] = tauVEffective;
            SUPPORT_SoA[BO_BATCH_SUPPORT_gInf_v * N + cellI] = vSteady;
            SUPPORT_SoA[BO_BATCH_SUPPORT_tau_w  * N + cellI] = tauWEffective;
            SUPPORT_SoA[BO_BATCH_SUPPORT_gInf_w * N + cellI] = wSteady;
            SUPPORT_SoA[BO_BATCH_SUPPORT_tau_s  * N + cellI] = tauSLocal;
            SUPPORT_SoA[BO_BATCH_SUPPORT_gInf_s * N + cellI] = sSteady;
            SUPPORT_SoA[BO_BATCH_SUPPORT_Iion   * N + cellI] = Jion;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
