/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    State-major SoA Batch evaluator for the TNNP 2004 ionic model.
    Architectural twin of `BuenoOrovio_2008Batch.H`.
\*---------------------------------------------------------------------------*/

#ifndef TNNP_2004Batch_H
#define TNNP_2004Batch_H

#include "TNNP_2004.H"
#include "TNNP_2004Names.H"
#include "StimulusProtocolPOD.H"

namespace Foam
{

// --- Compact support index map -----------------------------------------
//
// Order is intentionally kept compact and stable (RL gates first, Iion last)
// so the dispatch tables in `TNNPBatched.C::TNNPRushLarsenDispatch` resolve via
// constant-time indexing. Do not reorder without updating that table.
enum TNNP_BATCH_SUPPORT_INDEX
{
    TNNP_BATCH_SUPPORT_tau_xr1,
    TNNP_BATCH_SUPPORT_gInf_xr1,
    TNNP_BATCH_SUPPORT_tau_xr2,
    TNNP_BATCH_SUPPORT_gInf_xr2,
    TNNP_BATCH_SUPPORT_tau_xs,
    TNNP_BATCH_SUPPORT_gInf_xs,
    TNNP_BATCH_SUPPORT_tau_m,
    TNNP_BATCH_SUPPORT_gInf_m,
    TNNP_BATCH_SUPPORT_tau_h,
    TNNP_BATCH_SUPPORT_gInf_h,
    TNNP_BATCH_SUPPORT_tau_j,
    TNNP_BATCH_SUPPORT_gInf_j,
    TNNP_BATCH_SUPPORT_tau_d,
    TNNP_BATCH_SUPPORT_gInf_d,
    TNNP_BATCH_SUPPORT_tau_f,
    TNNP_BATCH_SUPPORT_gInf_f,
    TNNP_BATCH_SUPPORT_tau_s,
    TNNP_BATCH_SUPPORT_gInf_s,
    TNNP_BATCH_SUPPORT_tau_r,
    TNNP_BATCH_SUPPORT_gInf_r,
    TNNP_BATCH_SUPPORT_Iion_cm,

    NUM_TNNP_BATCH_SUPPORT
};


CARDIAC_HOST_DEVICE inline void TNNPwriteCompactSupportSoA
(
    const double* __restrict__ ALGEBRAIC,
    const int N,
    const int cellI,
    double* __restrict__ SUPPORT_SoA
)
{
    if (SUPPORT_SoA == nullptr)
    {
        return;
    }

    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr1  * N + cellI] = ALGEBRAIC[tau_xr1];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xr1 * N + cellI] = ALGEBRAIC[xr1_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr2  * N + cellI] = ALGEBRAIC[tau_xr2];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xr2 * N + cellI] = ALGEBRAIC[xr2_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xs   * N + cellI] = ALGEBRAIC[tau_xs];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xs  * N + cellI] = ALGEBRAIC[xs_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_m    * N + cellI] = ALGEBRAIC[tau_m];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_m   * N + cellI] = ALGEBRAIC[m_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_h    * N + cellI] = ALGEBRAIC[tau_h];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_h   * N + cellI] = ALGEBRAIC[h_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_j    * N + cellI] = ALGEBRAIC[tau_j];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_j   * N + cellI] = ALGEBRAIC[j_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_d    * N + cellI] = ALGEBRAIC[tau_d];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_d   * N + cellI] = ALGEBRAIC[d_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_f    * N + cellI] = ALGEBRAIC[tau_f];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_f   * N + cellI] = ALGEBRAIC[f_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_s    * N + cellI] = ALGEBRAIC[tau_s];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_s   * N + cellI] = ALGEBRAIC[s_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_r    * N + cellI] = ALGEBRAIC[tau_r];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_r   * N + cellI] = ALGEBRAIC[r_inf];
    SUPPORT_SoA[TNNP_BATCH_SUPPORT_Iion_cm  * N + cellI] = ALGEBRAIC[Iion_cm];
}


CARDIAC_HOST_DEVICE inline void TNNPwriteCompactSupport
(
    const double* __restrict__ ALGEBRAIC,
    double* __restrict__ SUPPORT
)
{
    if (SUPPORT == nullptr)
    {
        return;
    }

    SUPPORT[TNNP_BATCH_SUPPORT_tau_xr1]  = ALGEBRAIC[tau_xr1];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_xr1] = ALGEBRAIC[xr1_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_xr2]  = ALGEBRAIC[tau_xr2];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_xr2] = ALGEBRAIC[xr2_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_xs]   = ALGEBRAIC[tau_xs];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_xs]  = ALGEBRAIC[xs_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_m]    = ALGEBRAIC[tau_m];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_m]   = ALGEBRAIC[m_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_h]    = ALGEBRAIC[tau_h];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_h]   = ALGEBRAIC[h_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_j]    = ALGEBRAIC[tau_j];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_j]   = ALGEBRAIC[j_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_d]    = ALGEBRAIC[tau_d];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_d]   = ALGEBRAIC[d_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_f]    = ALGEBRAIC[tau_f];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_f]   = ALGEBRAIC[f_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_s]    = ALGEBRAIC[tau_s];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_s]   = ALGEBRAIC[s_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_r]    = ALGEBRAIC[tau_r];
    SUPPORT[TNNP_BATCH_SUPPORT_gInf_r]   = ALGEBRAIC[r_inf];
    SUPPORT[TNNP_BATCH_SUPPORT_Iion_cm]  = ALGEBRAIC[Iion_cm];
}


CARDIAC_HOST_DEVICE inline void TNNPComputeRatesBatchImpl
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ ALGEBRAIC_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVmWithinODESolver,
    const StimulusProtocolPOD& stimulus
)
{
    const int V_BASE = V * N;
    const int K_i_BASE = K_i * N;
    const int Na_i_BASE = Na_i * N;
    const int Ca_i_BASE = Ca_i * N;
    const int Xr1_BASE = Xr1 * N;
    const int Xr2_BASE = Xr2 * N;
    const int Xs_BASE = Xs * N;
    const int m_BASE = m * N;
    const int h_BASE = h * N;
    const int j_BASE = j * N;
    const int d_BASE = d * N;
    const int f_BASE = f * N;
    const int fCa_BASE = fCa * N;
    const int s_BASE = s * N;
    const int r_BASE = r * N;
    const int Ca_SR_BASE = Ca_SR * N;
    const int g_BASE = g * N;

    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double s_V = STATES_SoA[V_BASE + cellI];
        const double s_K_i = STATES_SoA[K_i_BASE + cellI];
        const double s_Na_i = STATES_SoA[Na_i_BASE + cellI];
        const double s_Ca_i = STATES_SoA[Ca_i_BASE + cellI];
        const double s_Xr1 = STATES_SoA[Xr1_BASE + cellI];
        const double s_Xr2 = STATES_SoA[Xr2_BASE + cellI];
        const double s_Xs = STATES_SoA[Xs_BASE + cellI];
        const double s_m = STATES_SoA[m_BASE + cellI];
        const double s_h = STATES_SoA[h_BASE + cellI];
        const double s_j = STATES_SoA[j_BASE + cellI];
        const double s_d = STATES_SoA[d_BASE + cellI];
        const double s_f = STATES_SoA[f_BASE + cellI];
        const double s_fCa = STATES_SoA[fCa_BASE + cellI];
        const double s_s = STATES_SoA[s_BASE + cellI];
        const double s_r = STATES_SoA[r_BASE + cellI];
        const double s_Ca_SR = STATES_SoA[Ca_SR_BASE + cellI];
        const double s_g = STATES_SoA[g_BASE + cellI];

        double STATES[NUM_STATES] = {0.0};
        double RATES[NUM_STATES] = {0.0};
        double ALGEBRAIC[NUM_ALGEBRAIC] = {0.0};

        STATES[V] = s_V;
        STATES[K_i] = s_K_i;
        STATES[Na_i] = s_Na_i;
        STATES[Ca_i] = s_Ca_i;
        STATES[Xr1] = s_Xr1;
        STATES[Xr2] = s_Xr2;
        STATES[Xs] = s_Xs;
        STATES[m] = s_m;
        STATES[h] = s_h;
        STATES[j] = s_j;
        STATES[d] = s_d;
        STATES[f] = s_f;
        STATES[fCa] = s_fCa;
        STATES[s] = s_s;
        STATES[r] = s_r;
        STATES[Ca_SR] = s_Ca_SR;
        STATES[g] = s_g;

        ALGEBRAIC[m_inf] = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
        ALGEBRAIC[alpha_m] = 1.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/5.00000));
        ALGEBRAIC[beta_m] = 0.100000/(1.00000+std::exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+std::exp((STATES[V] - 50.0000)/200.000));
        ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
        RATES[m] = (ALGEBRAIC[m_inf] - STATES[m])/ALGEBRAIC[tau_m];
        ALGEBRAIC[h_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
        ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*std::exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
        ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*std::exp( 0.0790000*STATES[V])+ 310000.*std::exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+std::exp((STATES[V]+10.6600)/- 11.1000))));
        ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
        RATES[h] = (ALGEBRAIC[h_inf] - STATES[h])/ALGEBRAIC[tau_h];
        ALGEBRAIC[j_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
        ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*STATES[V]) -  6.94800e-06*std::exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
        ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*STATES[V]))/(1.00000+std::exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*std::exp( 0.0570000*STATES[V]))/(1.00000+std::exp( - 0.100000*(STATES[V]+32.0000))));
        ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
        RATES[j] = (ALGEBRAIC[j_inf] - STATES[j])/ALGEBRAIC[tau_j];

        ALGEBRAIC[d_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/7.50000));
        ALGEBRAIC[alpha_d] = 1.40000/(1.00000+std::exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
        ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+std::exp((STATES[V]+5.00000)/5.00000));
        ALGEBRAIC[gamma_d] = 1.00000/(1.00000+std::exp((50.0000 - STATES[V])/20.0000));
        ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
        RATES[d] = (ALGEBRAIC[d_inf] - STATES[d])/ALGEBRAIC[tau_d];
        ALGEBRAIC[f_inf] = 1.00000/(1.00000+std::exp((STATES[V]+20.0000)/7.00000));
        ALGEBRAIC[tau_f] =  1125.00*std::exp(- std::pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - STATES[V])/10.0000));
        RATES[f] = (ALGEBRAIC[f_inf] - STATES[f])/ALGEBRAIC[tau_f];

        ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000325000, 8.00000));
        ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+std::exp((STATES[Ca_i] - 0.000500000)/0.000100000));
        ALGEBRAIC[beta_d] = 0.200000/(1.00000+std::exp((STATES[Ca_i] - 0.000750000)/0.000800000));
        ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
        ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
        RATES[fCa] = (ALGEBRAIC[fCa_inf]>STATES[fCa]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_fCa]);

        ALGEBRAIC[s_inf] = CONSTANTS[sInfPrefactor]
            / (1.00000 + std::exp((STATES[V] + CONSTANTS[sInfShift]) / CONSTANTS[sInfScale]));
        ALGEBRAIC[tau_s] =
              CONSTANTS[tauSGaussAmp]
              * std::exp(- std::pow(STATES[V] + CONSTANTS[tauSGaussShift], 2.00000)
                    / CONSTANTS[tauSGaussWidth])
            + CONSTANTS[tauSSigmoidAmp]
              / (1.00000 + std::exp((STATES[V] + CONSTANTS[tauSSigmoidShift])
                               / CONSTANTS[tauSSigmoidScale]))
            + CONSTANTS[tauSOffset];
        RATES[s] = (ALGEBRAIC[s_inf] - STATES[s])/ALGEBRAIC[tau_s];
        ALGEBRAIC[r_inf] = 1.00000/(1.00000+std::exp((20.0000 - STATES[V])/6.00000));
        ALGEBRAIC[tau_r] =  9.50000*std::exp(- std::pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
        RATES[r] = (ALGEBRAIC[r_inf] - STATES[r])/ALGEBRAIC[tau_r];

        ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+std::exp((- 26.0000 - STATES[V])/7.00000));
        ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+std::exp((- 45.0000 - STATES[V])/10.0000));
        ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+std::exp((STATES[V]+30.0000)/11.5000));
        ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
        RATES[Xr1] = (ALGEBRAIC[xr1_inf] - STATES[Xr1])/ALGEBRAIC[tau_xr1];
        ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+std::exp((STATES[V]+88.0000)/24.0000));
        ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/20.0000));
        ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
        ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
        RATES[Xr2] = (ALGEBRAIC[xr2_inf] - STATES[Xr2])/ALGEBRAIC[tau_xr2];

        ALGEBRAIC[xs_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/14.0000));
        ALGEBRAIC[alpha_xs] = 1100.00/ std::pow((1.00000+std::exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
        ALGEBRAIC[beta_xs] = 1.00000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
        ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
        RATES[Xs] = (ALGEBRAIC[xs_inf] - STATES[Xs])/ALGEBRAIC[tau_xs];

        ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 16.0000)));
        ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
        RATES[g] = (ALGEBRAIC[g_inf]>STATES[g]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_g]);


        ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*std::exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T]))+ 0.0353000*std::exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])));
        ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])/CONSTANTS[F])*std::log(CONSTANTS[Na_o]/STATES[Na_i]);
        ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*std::pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
        ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
        ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( std::exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T]))*std::pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T]))*std::pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (std::pow(CONSTANTS[Km_Nai], 3.00000)+std::pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T]))));
        RATES[Na_i] = ( - (ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ 3.00000*ALGEBRAIC[i_NaK]+ 3.00000*ALGEBRAIC[i_NaCa])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);

        ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])/CONSTANTS[F])*std::log(CONSTANTS[K_o]/STATES[K_i]);
        ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+std::exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
        ALGEBRAIC[beta_K1] = ( 3.00000*std::exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*std::exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+std::exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
        ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
        ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);

        ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);

        ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
        ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])/CONSTANTS[F])*std::log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
        ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*std::pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
        ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*std::pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T]))*( STATES[Ca_i]*std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])) -  0.341000*CONSTANTS[Ca_o]))/(std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])) - 1.00000);
        ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[CONSTANTS_INDEX::T])/CONSTANTS[F])*std::log(CONSTANTS[Ca_o]/STATES[Ca_i]);
        ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
        ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+std::exp((25.0000 - STATES[V])/5.98000));
        ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

        ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                        ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                        ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                        ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

        ALGEBRAIC[Istim] = computeStimulusPOD(t, stimulus);

        RATES[V] = 0.0;
        if (solveVmWithinODESolver)
        {
            RATES[V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
        }

        RATES[K_i] = ( - ((ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_p_K]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[i_NaK])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);
        ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*std::pow(STATES[Ca_SR], 2.00000))/(std::pow(CONSTANTS[b_rel], 2.00000)+std::pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
        ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+std::pow(CONSTANTS[K_up], 2.00000)/std::pow(STATES[Ca_i], 2.00000));
        ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
        ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
        ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/std::pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
        RATES[Ca_i] =  ALGEBRAIC[ddt_Ca_i_total]*ALGEBRAIC[f_JCa_i_free];
        ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
        ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/std::pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));
        RATES[Ca_SR] =  ALGEBRAIC[ddt_Ca_sr_total]*ALGEBRAIC[f_JCa_sr_free];

        RATES_SoA[V_BASE + cellI] = RATES[V];
        RATES_SoA[K_i_BASE + cellI] = RATES[K_i];
        RATES_SoA[Na_i_BASE + cellI] = RATES[Na_i];
        RATES_SoA[Ca_i_BASE + cellI] = RATES[Ca_i];
        RATES_SoA[Xr1_BASE + cellI] = RATES[Xr1];
        RATES_SoA[Xr2_BASE + cellI] = RATES[Xr2];
        RATES_SoA[Xs_BASE + cellI] = RATES[Xs];
        RATES_SoA[m_BASE + cellI] = RATES[m];
        RATES_SoA[h_BASE + cellI] = RATES[h];
        RATES_SoA[j_BASE + cellI] = RATES[j];
        RATES_SoA[d_BASE + cellI] = RATES[d];
        RATES_SoA[f_BASE + cellI] = RATES[f];
        RATES_SoA[fCa_BASE + cellI] = RATES[fCa];
        RATES_SoA[s_BASE + cellI] = RATES[s];
        RATES_SoA[r_BASE + cellI] = RATES[r];
        RATES_SoA[Ca_SR_BASE + cellI] = RATES[Ca_SR];
        RATES_SoA[g_BASE + cellI] = RATES[g];

        if (ALGEBRAIC_SoA != nullptr)
        {
            for (int algI = 0; algI < NUM_ALGEBRAIC; ++algI)
            {
                ALGEBRAIC_SoA[algI*N + cellI] = ALGEBRAIC[algI];
            }
        }

        TNNPwriteCompactSupportSoA(ALGEBRAIC, N, cellI, SUPPORT_SoA);
    }
}


CARDIAC_HOST_DEVICE inline void TNNPComputeRatesBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVmWithinODESolver,
    const StimulusProtocolPOD& stimulus
)
{
    TNNPComputeRatesBatchImpl
    (
        t,
        CONSTANTS,
        N,
        beginCell,
        endCell,
        STATES_SoA,
        RATES_SoA,
        nullptr,
        SUPPORT_SoA,
        solveVmWithinODESolver,
        stimulus
    );
}


CARDIAC_HOST_DEVICE inline void TNNPComputeRatesFullBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ ALGEBRAIC_SoA,
    const bool solveVmWithinODESolver,
    const StimulusProtocolPOD& stimulus
)
{
    TNNPComputeRatesBatchImpl
    (
        t,
        CONSTANTS,
        N,
        beginCell,
        endCell,
        STATES_SoA,
        RATES_SoA,
        ALGEBRAIC_SoA,
        nullptr,
        solveVmWithinODESolver,
        stimulus
    );
}


// --- CPU-only compact-support bridge -----------------------------------
//
// TNNPcomputeRatesCompact is a *host-only* delegating shim.
// It accepts a `Foam::StimulusProtocol&` argument (a Foam type) and
// therefore cannot be compiled for or called from device code.
//
// Purpose: bridges `evaluateHotPathState()` (the per-cell virtual-dispatch
// path used by Rush-Larsen and utilities mode) to the scalar `TNNPcomputeRates`
// evaluator, projecting its ALGEBRAIC output into the compact SUPPORT layout.
//
// The SoA batch path (`solveBatched` / CUDA kernel) calls
// `TNNPComputeRatesBatch` directly and never reaches this function.
inline void
TNNPcomputeRatesCompact
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* SUPPORT,
    bool solveVmWithinODESolver,
    const Foam::StimulusProtocol& stimulus
)
{
    double ALGEBRAIC[NUM_ALGEBRAIC] = {0.0};

    TNNPcomputeRates
    (
        VOI,
        CONSTANTS,
        RATES,
        STATES,
        ALGEBRAIC,
        solveVmWithinODESolver,
        stimulus
    );

    TNNPwriteCompactSupport(ALGEBRAIC, SUPPORT);
}

} // End namespace Foam

#endif
