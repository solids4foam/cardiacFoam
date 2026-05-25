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
    TNNP_BATCH_SUPPORT_tau_xr2,
    TNNP_BATCH_SUPPORT_tau_xs,
    TNNP_BATCH_SUPPORT_tau_m,
    TNNP_BATCH_SUPPORT_tau_h,
    TNNP_BATCH_SUPPORT_tau_j,
    TNNP_BATCH_SUPPORT_tau_d,
    TNNP_BATCH_SUPPORT_tau_f,
    TNNP_BATCH_SUPPORT_tau_s,
    TNNP_BATCH_SUPPORT_tau_r,
    TNNP_BATCH_SUPPORT_Iion_cm,

    NUM_TNNP_BATCH_SUPPORT
};


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
    // Stimulus is loop-invariant - compute once.
    const double Istim = computeStimulusPOD(t, stimulus);

    // Hoist constant reads outside the cell loop.
    const double c_R = CONSTANTS[R];
    const double c_T = CONSTANTS[CONSTANTS_INDEX::T];
    const double c_F = CONSTANTS[F];
    const double c_Cm = CONSTANTS[Cm];
    const double c_V_c = CONSTANTS[V_c];
    const double c_P_kna = CONSTANTS[P_kna];
    const double c_K_o = CONSTANTS[K_o];
    const double c_Na_o = CONSTANTS[Na_o];
    const double c_Ca_o = CONSTANTS[Ca_o];
    const double c_g_K1 = CONSTANTS[g_K1];
    const double c_g_Kr = CONSTANTS[g_Kr];
    const double c_g_Ks = CONSTANTS[g_Ks];
    const double c_g_Na = CONSTANTS[g_Na];
    const double c_g_bNa = CONSTANTS[g_bNa];
    const double c_g_CaL = CONSTANTS[g_CaL];
    const double c_g_bCa = CONSTANTS[g_bCa];
    const double c_g_to = CONSTANTS[g_to];
    const double c_P_NaK = CONSTANTS[P_NaK];
    const double c_K_mk = CONSTANTS[K_mk];
    const double c_K_mNa = CONSTANTS[K_mNa];
    const double c_K_NaCa = CONSTANTS[K_NaCa];
    const double c_K_sat = CONSTANTS[K_sat];
    const double c_alpha_NCX = CONSTANTS[alpha_NCX];
    const double c_gamma_NCX = CONSTANTS[gamma_NCX];
    const double c_Km_Ca = CONSTANTS[Km_Ca];
    const double c_Km_Nai = CONSTANTS[Km_Nai];
    const double c_g_pCa = CONSTANTS[g_pCa];
    const double c_K_pCa = CONSTANTS[K_pCa];
    const double c_g_pK = CONSTANTS[g_pK];
    const double c_tau_g = CONSTANTS[tau_g];
    const double c_a_rel = CONSTANTS[a_rel];
    const double c_b_rel = CONSTANTS[b_rel];
    const double c_c_rel = CONSTANTS[c_rel];
    const double c_K_up = CONSTANTS[K_up];
    const double c_V_leak = CONSTANTS[V_leak];
    const double c_Vmax_up = CONSTANTS[Vmax_up];
    const double c_Buf_c = CONSTANTS[Buf_c];
    const double c_K_buf_c = CONSTANTS[K_buf_c];
    const double c_Buf_sr = CONSTANTS[Buf_sr];
    const double c_K_buf_sr = CONSTANTS[K_buf_sr];
    const double c_V_sr = CONSTANTS[V_sr];
    const double c_tau_fCa = CONSTANTS[tau_fCa];

    // SoA strides for each state.
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

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        // Load all states.
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

        // INa
        const double m_inf = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - s_V)/9.03000), 2.00000);
        const double alpha_m = 1.00000/(1.00000+std::exp((- 60.0000 - s_V)/5.00000));
        const double beta_m = 0.100000/(1.00000+std::exp((s_V+35.0000)/5.00000))+0.100000/(1.00000+std::exp((s_V - 50.0000)/200.000));
        const double tau_m =  1.00000*alpha_m*beta_m;
        const double r_m = (m_inf - s_m)/tau_m;

        const double h_inf = 1.00000/std::pow(1.00000+std::exp((s_V+71.5500)/7.43000), 2.00000);
        const double alpha_h = (s_V<- 40.0000 ?  0.0570000*std::exp(- (s_V+80.0000)/6.80000) : 0.00000);
        const double beta_h = (s_V<- 40.0000 ?  2.70000*std::exp( 0.0790000*s_V)+ 310000.*std::exp( 0.348500*s_V) : 0.770000/( 0.130000*(1.00000+std::exp((s_V+10.6600)/- 11.1000))));
        const double tau_h = 1.00000/(alpha_h+beta_h);
        const double r_h = (h_inf - s_h)/tau_h;

        const double j_inf = 1.00000/std::pow(1.00000+std::exp((s_V+71.5500)/7.43000), 2.00000);
        const double alpha_j = (s_V<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*s_V) -  6.94800e-06*std::exp( - 0.0439100*s_V))*(s_V+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(s_V+79.2300))) : 0.00000);
        const double beta_j = (s_V<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*s_V))/(1.00000+std::exp( - 0.137800*(s_V+40.1400))) : ( 0.600000*std::exp( 0.0570000*s_V))/(1.00000+std::exp( - 0.100000*(s_V+32.0000))));
        const double tau_j = 1.00000/(alpha_j+beta_j);
        const double r_j = (j_inf - s_j)/tau_j;

        // ICaL
        const double d_inf = 1.00000/(1.00000+std::exp((- 5.00000 - s_V)/7.50000));
        const double alpha_d = 1.40000/(1.00000+std::exp((- 35.0000 - s_V)/13.0000))+0.250000;
        const double gama_fCa = 1.40000/(1.00000+std::exp((s_V+5.00000)/5.00000));
        const double gamma_d = 1.00000/(1.00000+std::exp((50.0000 - s_V)/20.0000));
        const double tau_d =  1.00000*alpha_d*gama_fCa+gamma_d;
        const double r_d = (d_inf - s_d)/tau_d;

        const double f_inf = 1.00000/(1.00000+std::exp((s_V+20.0000)/7.00000));
        const double tau_f =  1125.00*std::exp(- std::pow(s_V+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - s_V)/10.0000));
        const double r_f = (f_inf - s_f)/tau_f;

        const double alpha_fCa = 1.00000/(1.00000+std::pow(s_Ca_i/0.000325000, 8.00000));
        const double beta_fCa = 0.100000/(1.00000+std::exp((s_Ca_i - 0.000500000)/0.000100000));
        const double beta_d = 0.200000/(1.00000+std::exp((s_Ca_i - 0.000750000)/0.000800000));
        const double fCa_inf = (alpha_fCa+beta_fCa+beta_d+0.230000)/1.46000;
        const double d_fCa_val = (fCa_inf - s_fCa)/c_tau_fCa;
        const double r_fCa = (fCa_inf>s_fCa&&s_V>- 60.0000 ? 0.00000 : d_fCa_val);

        // Ito
        const double s_inf = CONSTANTS[sInfPrefactor]
            / (1.00000 + std::exp((s_V + CONSTANTS[sInfShift]) / CONSTANTS[sInfScale]));
        const double tau_s =
              CONSTANTS[tauSGaussAmp]
              * std::exp(- std::pow(s_V + CONSTANTS[tauSGaussShift], 2.00000)
                    / CONSTANTS[tauSGaussWidth])
            + CONSTANTS[tauSSigmoidAmp]
              / (1.00000 + std::exp((s_V + CONSTANTS[tauSSigmoidShift])
                               / CONSTANTS[tauSSigmoidScale]))
            + CONSTANTS[tauSOffset];
        const double r_s = (s_inf - s_s)/tau_s;

        const double r_inf = 1.00000/(1.00000+std::exp((20.0000 - s_V)/6.00000));
        const double tau_r =  9.50000*std::exp(- std::pow(s_V+40.0000, 2.00000)/1800.00)+0.800000;
        const double r_r = (r_inf - s_r)/tau_r;

        // IKr
        const double xr1_inf = 1.00000/(1.00000+std::exp((- 26.0000 - s_V)/7.00000));
        const double alpha_xr1 = 450.000/(1.00000+std::exp((- 45.0000 - s_V)/10.0000));
        const double beta_xr1 = 6.00000/(1.00000+std::exp((s_V+30.0000)/11.5000));
        const double tau_xr1 =  1.00000*alpha_xr1*beta_xr1;
        const double r_Xr1 = (xr1_inf - s_Xr1)/tau_xr1;

        const double xr2_inf = 1.00000/(1.00000+std::exp((s_V+88.0000)/24.0000));
        const double alpha_xr2 = 3.00000/(1.00000+std::exp((- 60.0000 - s_V)/20.0000));
        const double beta_xr2 = 1.12000/(1.00000+std::exp((s_V - 60.0000)/20.0000));
        const double tau_xr2 =  1.00000*alpha_xr2*beta_xr2;
        const double r_Xr2 = (xr2_inf - s_Xr2)/tau_xr2;

        // IKs
        const double xs_inf = 1.00000/(1.00000+std::exp((- 5.00000 - s_V)/14.0000));
        const double alpha_xs = 1100.00/std::sqrt(1.00000+std::exp((- 10.0000 - s_V)/6.00000));
        const double beta_xs = 1.00000/(1.00000+std::exp((s_V - 60.0000)/20.0000));
        const double tau_xs =  1.00000*alpha_xs*beta_xs;
        const double r_Xs = (xs_inf - s_Xs)/tau_xs;

        // SR calcium RyR gate g
        const double g_inf = (s_Ca_i<0.000350000 ? 1.00000/(1.00000+std::pow(s_Ca_i/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(s_Ca_i/0.000350000, 16.0000)));
        const double d_g_val = (g_inf - s_g)/c_tau_g;
        const double r_g = (g_inf>s_g&&s_V>- 60.0000 ? 0.00000 : d_g_val);

        // Other currents
        const double i_NaK = (( (( c_P_NaK*c_K_o)/(c_K_o+c_K_mk))*s_Na_i)/(s_Na_i+c_K_mNa))/(1.00000+ 0.124500*std::exp(( - 0.100000*s_V*c_F)/( c_R*c_T))+ 0.0353000*std::exp(( - s_V*c_F)/( c_R*c_T)));
        const double E_Na =  (( c_R*c_T)/c_F)*std::log(c_Na_o/s_Na_i);
        const double i_Na =  c_g_Na*std::pow(s_m, 3.00000)*s_h*s_j*(s_V - E_Na);
        const double i_b_Na =  c_g_bNa*(s_V - E_Na);
        const double i_NaCa = ( c_K_NaCa*( std::exp(( c_gamma_NCX*s_V*c_F)/( c_R*c_T))*std::pow(s_Na_i, 3.00000)*c_Ca_o -  std::exp(( (c_gamma_NCX - 1.00000)*s_V*c_F)/( c_R*c_T))*std::pow(c_Na_o, 3.00000)*s_Ca_i*c_alpha_NCX))/( (std::pow(c_Km_Nai, 3.00000)+std::pow(c_Na_o, 3.00000))*(c_Km_Ca+c_Ca_o)*(1.00000+ c_K_sat*std::exp(( (c_gamma_NCX - 1.00000)*s_V*c_F)/( c_R*c_T))));
        const double r_Na_i = ( - (i_Na+i_b_Na+ 3.00000*i_NaK+ 3.00000*i_NaCa)*c_Cm)/( c_V_c*c_F);

        const double E_K =  (( c_R*c_T)/c_F)*std::log(c_K_o/s_K_i);
        const double alpha_K1 = 0.100000/(1.00000+std::exp( 0.0600000*((s_V - E_K) - 200.000)));
        const double beta_K1 = ( 3.00000*std::exp( 0.000200000*((s_V - E_K)+100.000))+ 1.00000*std::exp( 0.100000*((s_V - E_K) - 10.0000)))/(1.00000+std::exp( - 0.500000*(s_V - E_K)));
        const double xK1_inf = alpha_K1/(alpha_K1+beta_K1);
        const double i_K1 =  c_g_K1*xK1_inf* std::sqrt(c_K_o/5.40000)*(s_V - E_K);

        const double i_to =  c_g_to*s_r*s_s*(s_V - E_K);

        const double i_Kr =  c_g_Kr* std::sqrt(c_K_o/5.40000)*s_Xr1*s_Xr2*(s_V - E_K);
        const double E_Ks =  (( c_R*c_T)/c_F)*std::log((c_K_o+ c_P_kna*c_Na_o)/(s_K_i+ c_P_kna*s_Na_i));
        const double i_Ks =  c_g_Ks*std::pow(s_Xs, 2.00000)*(s_V - E_Ks);

        const double exp2VFRT = std::exp((2.00000*s_V*c_F)/(c_R*c_T));
        const double i_CaL = (((c_g_CaL*s_d*s_f*s_fCa*4.00000*s_V*std::pow(c_F, 2.00000))/(c_R*c_T)) * (s_Ca_i*exp2VFRT - 0.341000*c_Ca_o)) / (exp2VFRT - 1.00000);

        const double E_Ca =  (( 0.500000*c_R*c_T)/c_F)*std::log(c_Ca_o/s_Ca_i);
        const double i_b_Ca =  c_g_bCa*(s_V - E_Ca);
        const double i_p_K = ( c_g_pK*(s_V - E_K))/(1.00000+std::exp((25.0000 - s_V)/5.98000));
        const double i_p_Ca = ( c_g_pCa*s_Ca_i)/(s_Ca_i+c_K_pCa);

        const double Iion_cm_val = i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_b_Na+i_NaCa+i_b_Ca+i_p_K+i_p_Ca;

        const double r_V = solveVmWithinODESolver ? (-Iion_cm_val - Istim) : 0.0;

        const double r_K_i = ( - ((i_K1+i_to+i_Kr+i_Ks+i_p_K+Istim) -  2.00000*i_NaK)*c_Cm)/( c_V_c*c_F);
        const double i_rel =  (( c_a_rel*std::pow(s_Ca_SR, 2.00000))/(std::pow(c_b_rel, 2.00000)+std::pow(s_Ca_SR, 2.00000))+c_c_rel)*s_d*s_g;
        const double i_up = c_Vmax_up/(1.00000+std::pow(c_K_up, 2.00000)/std::pow(s_Ca_i, 2.00000));
        const double i_leak =  c_V_leak*(s_Ca_SR - s_Ca_i);
        const double ddt_Ca_i_total = (( (- ((i_CaL+i_b_Ca+i_p_Ca) -  2.00000*i_NaCa)/( 2.00000*c_V_c*c_F))*c_Cm+i_leak) - i_up)+i_rel;
        const double f_JCa_i_free = 1.00000/(1.00000+( c_Buf_c*c_K_buf_c)/std::pow(s_Ca_i+c_K_buf_c, 2.00000));
        const double r_Ca_i =  ddt_Ca_i_total*f_JCa_i_free;
        const double ddt_Ca_sr_total =  (c_V_c/c_V_sr)*(i_up - (i_rel+i_leak));
        const double f_JCa_sr_free = 1.00000/(1.00000+( c_Buf_sr*c_K_buf_sr)/std::pow(s_Ca_SR+c_K_buf_sr, 2.00000));
        const double r_Ca_SR =  ddt_Ca_sr_total*f_JCa_sr_free;

        // Write rates back to SoA
        RATES_SoA[V_BASE + cellI] = r_V;
        RATES_SoA[K_i_BASE + cellI] = r_K_i;
        RATES_SoA[Na_i_BASE + cellI] = r_Na_i;
        RATES_SoA[Ca_i_BASE + cellI] = r_Ca_i;
        RATES_SoA[Xr1_BASE + cellI] = r_Xr1;
        RATES_SoA[Xr2_BASE + cellI] = r_Xr2;
        RATES_SoA[Xs_BASE + cellI] = r_Xs;
        RATES_SoA[m_BASE + cellI] = r_m;
        RATES_SoA[h_BASE + cellI] = r_h;
        RATES_SoA[j_BASE + cellI] = r_j;
        RATES_SoA[d_BASE + cellI] = r_d;
        RATES_SoA[f_BASE + cellI] = r_f;
        RATES_SoA[fCa_BASE + cellI] = r_fCa;
        RATES_SoA[s_BASE + cellI] = r_s;
        RATES_SoA[r_BASE + cellI] = r_r;
        RATES_SoA[Ca_SR_BASE + cellI] = r_Ca_SR;
        RATES_SoA[g_BASE + cellI] = r_g;

        // Write optional compact-support back to SoA
        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr1 * N + cellI] = tau_xr1;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr2 * N + cellI] = tau_xr2;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xs  * N + cellI] = tau_xs;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_m   * N + cellI] = tau_m;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_h   * N + cellI] = tau_h;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_j   * N + cellI] = tau_j;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_d   * N + cellI] = tau_d;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_f   * N + cellI] = tau_f;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_s   * N + cellI] = tau_s;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_r   * N + cellI] = tau_r;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_Iion_cm * N + cellI] = Iion_cm_val;
        }
    }
}


// --- Compact evaluator -------------------------------------------------
//
// Per-cell delegating shim: runs the full TNNPcomputeRates and projects the
// result into the compact SUPPORT layout. TNNP carries this shim (Bueno/Gaur
// do not) because TNNP_2004.H exposes two per-cell evaluators
// (TNNPcomputeRates and TNNPcomputeVariables) with different surface shapes
// from what the SoA batched loop produces. Future work may replace this
// delegating body with a directly-factored compact expression set.
//
// The `tissueFlag` argument is accepted for call-site compatibility with
// the other per-cell evaluators but is unused here (tissue heterogeneity
// is baked into CONSTANTS by TNNPinitConsts).

// --- CPU-only bridging shim -------------------------------------------
//
// TNNPcomputeRatesCompact is a *host-only* delegating shim.
// It accepts a `Foam::StimulusProtocol&` argument (a Foam type) and
// therefore cannot be compiled for or called from device code.
//
// Purpose: bridges `evaluateHotPathState()` (the per-cell virtual-dispatch
// fallback path used by Rush-Larsen and the utilities mode) to the scalar
// `TNNPcomputeRates` evaluator, projecting its ALGEBRAIC output into the
// compact SUPPORT layout expected by the hot-path integrator.
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

    SUPPORT[TNNP_BATCH_SUPPORT_tau_xr1] = ALGEBRAIC[tau_xr1];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_xr2] = ALGEBRAIC[tau_xr2];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_xs]  = ALGEBRAIC[tau_xs];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_m]   = ALGEBRAIC[tau_m];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_h]   = ALGEBRAIC[tau_h];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_j]   = ALGEBRAIC[tau_j];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_d]   = ALGEBRAIC[tau_d];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_f]   = ALGEBRAIC[tau_f];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_s]   = ALGEBRAIC[tau_s];
    SUPPORT[TNNP_BATCH_SUPPORT_tau_r]   = ALGEBRAIC[tau_r];
    SUPPORT[TNNP_BATCH_SUPPORT_Iion_cm] = ALGEBRAIC[Iion_cm];
}

} // End namespace Foam

#endif
