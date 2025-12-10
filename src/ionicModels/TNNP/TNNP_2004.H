/*
   There are a total of 69 entries in the algebraic variable array.
   There are a total of 17 entries in each of the rate and state variable arrays.
   There are a total of 46 entries in the constant variable array +1 for the flag that defines the tissue being modeled.

 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_millimole).
 * CONSTANTS[3] is Cm in component membrane (picoF).
 * CONSTANTS[4] is V_c in component membrane (micrometre3).
 * ALGEBRAIC[50] is i_K1 in component inward_rectifier_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[57] is i_to in component transient_outward_current (picoA_per_picoF).
 * ALGEBRAIC[51] is i_Kr in component rapid_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[52] is i_Ks in component slow_time_dependent_potassium_current (picoA_per_picoF).
 * ALGEBRAIC[55] is i_CaL in component L_type_Ca_current (picoA_per_picoF).
 * ALGEBRAIC[58] is i_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[53] is i_Na in component fast_sodium_current (picoA_per_picoF).
 * ALGEBRAIC[54] is i_b_Na in component sodium_background_current (picoA_per_picoF).
 * ALGEBRAIC[59] is i_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * ALGEBRAIC[56] is i_b_Ca in component calcium_background_current (picoA_per_picoF).
 * ALGEBRAIC[61] is i_p_K in component potassium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[60] is i_p_Ca in component calcium_pump_current (picoA_per_picoF).
 * ALGEBRAIC[0] is i_Stim in component membrane (picoA_per_picoF).
 * CONSTANTS[5] is stim_start in component membrane (millisecond).
 * CONSTANTS[6] is stim_period in component membrane (millisecond).
 * CONSTANTS[7] is stim_duration in component membrane (millisecond).
 * CONSTANTS[8] is stim_amplitude in component membrane (picoA_per_picoF).
 * ALGEBRAIC[13] is E_Na in component reversal_potentials (millivolt).
 * ALGEBRAIC[26] is E_K in component reversal_potentials (millivolt).
 * ALGEBRAIC[35] is E_Ks in component reversal_potentials (millivolt).
 * ALGEBRAIC[44] is E_Ca in component reversal_potentials (millivolt).
 * CONSTANTS[9] is P_kna in component reversal_potentials (dimensionless).
 * CONSTANTS[10] is K_o in component potassium_dynamics (millimolar).
 * CONSTANTS[11] is Na_o in component sodium_dynamics (millimolar).
 * STATES[1] is K_i in component potassium_dynamics (millimolar).
 * STATES[2] is Na_i in component sodium_dynamics (millimolar).
 * CONSTANTS[12] is Ca_o in component calcium_dynamics (millimolar).
 * STATES[3] is Ca_i in component calcium_dynamics (millimolar).
 * CONSTANTS[13] is g_K1 in component inward_rectifier_potassium_current (nanoS_per_picoF).
 * ALGEBRAIC[49] is xK1_inf in component inward_rectifier_potassium_current (dimensionless).
 * ALGEBRAIC[47] is alpha_K1 in component inward_rectifier_potassium_current (per_millisecond).
 * ALGEBRAIC[48] is beta_K1 in component inward_rectifier_potassium_current (per_millisecond).
 * CONSTANTS[14] is g_Kr in component rapid_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[4] is Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * STATES[5] is Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[1] is xr1_inf in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[27] is beta_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * ALGEBRAIC[36] is tau_xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (millisecond).
 * ALGEBRAIC[2] is xr2_inf in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[15] is alpha_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[28] is beta_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * ALGEBRAIC[37] is tau_xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (millisecond).
 * CONSTANTS[15] is g_Ks in component slow_time_dependent_potassium_current (nanoS_per_picoF).
 * STATES[6] is Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[3] is xs_inf in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[29] is beta_xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * ALGEBRAIC[38] is tau_xs in component slow_time_dependent_potassium_current_Xs_gate (millisecond).
 * CONSTANTS[16] is g_Na in component fast_sodium_current (nanoS_per_picoF).
 * STATES[7] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[8] is h in component fast_sodium_current_h_gate (dimensionless).
 * STATES[9] is j in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[4] is m_inf in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[17] is alpha_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[30] is beta_m in component fast_sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[39] is tau_m in component fast_sodium_current_m_gate (millisecond).
 * ALGEBRAIC[5] is h_inf in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[18] is alpha_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[31] is beta_h in component fast_sodium_current_h_gate (per_millisecond).
 * ALGEBRAIC[40] is tau_h in component fast_sodium_current_h_gate (millisecond).
 * ALGEBRAIC[6] is j_inf in component fast_sodium_current_j_gate (dimensionless).
 * ALGEBRAIC[19] is alpha_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[32] is beta_j in component fast_sodium_current_j_gate (per_millisecond).
 * ALGEBRAIC[41] is tau_j in component fast_sodium_current_j_gate (millisecond).
 * CONSTANTS[17] is g_bna in component sodium_background_current (nanoS_per_picoF).
 * CONSTANTS[18] is g_CaL in component L_type_Ca_current (litre_per_farad_second).
 * STATES[10] is d in component L_type_Ca_current_d_gate (dimensionless).
 * STATES[11] is f in component L_type_Ca_current_f_gate (dimensionless).
 * STATES[12] is fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[7] is d_inf in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[20] is alpha_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[33] is beta_d in component L_type_Ca_current_d_gate (dimensionless).
 * ALGEBRAIC[42] is gamma_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[45] is tau_d in component L_type_Ca_current_d_gate (millisecond).
 * ALGEBRAIC[8] is f_inf in component L_type_Ca_current_f_gate (dimensionless).
 * ALGEBRAIC[21] is tau_f in component L_type_Ca_current_f_gate (millisecond).
 * ALGEBRAIC[9] is alpha_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[22] is beta_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[34] is gama_fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * ALGEBRAIC[43] is fCa_inf in component L_type_Ca_current_fCa_gate (dimensionless).
 * CONSTANTS[45] is tau_fCa in component L_type_Ca_current_fCa_gate (millisecond).
 * ALGEBRAIC[46] is d_fCa in component L_type_Ca_current_fCa_gate (per_millisecond).
 * CONSTANTS[19] is g_bca in component calcium_background_current (nanoS_per_picoF).
 * CONSTANTS[20] is g_to in component transient_outward_current (nanoS_per_picoF).
 * STATES[13] is s in component transient_outward_current_s_gate (dimensionless).
 * STATES[14] is r in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[10] is s_inf in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[23] is tau_s in component transient_outward_current_s_gate (millisecond).
 * ALGEBRAIC[11] is r_inf in component transient_outward_current_r_gate (dimensionless).
 * ALGEBRAIC[24] is tau_r in component transient_outward_current_r_gate (millisecond).
 * CONSTANTS[21] is P_NaK in component sodium_potassium_pump_current (picoA_per_picoF).
 * CONSTANTS[22] is K_mk in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[23] is K_mNa in component sodium_potassium_pump_current (millimolar).
 * CONSTANTS[24] is K_NaCa in component sodium_calcium_exchanger_current (picoA_per_picoF).
 * CONSTANTS[25] is K_sat in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[26] is alpha in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[27] is gamma in component sodium_calcium_exchanger_current (dimensionless).
 * CONSTANTS[28] is Km_Ca in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[29] is Km_Nai in component sodium_calcium_exchanger_current (millimolar).
 * CONSTANTS[30] is g_pCa in component calcium_pump_current (picoA_per_picoF).
 * CONSTANTS[31] is K_pCa in component calcium_pump_current (millimolar).
 * CONSTANTS[32] is g_pK in component potassium_pump_current (nanoS_per_picoF).
 * STATES[15] is Ca_SR in component calcium_dynamics (millimolar).
 * ALGEBRAIC[62] is i_rel in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[63] is i_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[64] is i_leak in component calcium_dynamics (millimolar_per_millisecond).
 * STATES[16] is g in component calcium_dynamics (dimensionless).
 * CONSTANTS[33] is tau_g in component calcium_dynamics (millisecond).
 * ALGEBRAIC[12] is g_inf in component calcium_dynamics (dimensionless).
 * CONSTANTS[34] is a_rel in component calcium_dynamics (millimolar_per_millisecond).
 * CONSTANTS[35] is b_rel in component calcium_dynamics (millimolar).
 * CONSTANTS[36] is c_rel in component calcium_dynamics (millimolar_per_millisecond).
 * CONSTANTS[37] is K_up in component calcium_dynamics (millimolar).
 * CONSTANTS[38] is V_leak in component calcium_dynamics (per_millisecond).
 * CONSTANTS[39] is Vmax_up in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[65] is ddt_Ca_i_total in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[66] is ddt_Ca_sr_total in component calcium_dynamics (millimolar_per_millisecond).
 * ALGEBRAIC[67] is f_JCa_i_free in component calcium_dynamics (dimensionless).
 * ALGEBRAIC[68] is f_JCa_sr_free in component calcium_dynamics (dimensionless).
 * CONSTANTS[40] is Buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[41] is K_buf_c in component calcium_dynamics (millimolar).
 * CONSTANTS[42] is Buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[43] is K_buf_sr in component calcium_dynamics (millimolar).
 * CONSTANTS[44] is V_sr in component calcium_dynamics (micrometre3).
 * ALGEBRAIC[25] is d_g in component calcium_dynamics (per_millisecond).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt Xr1 in component rapid_time_dependent_potassium_current_Xr1_gate (dimensionless).
 * RATES[5] is d/dt Xr2 in component rapid_time_dependent_potassium_current_Xr2_gate (dimensionless).
 * RATES[6] is d/dt Xs in component slow_time_dependent_potassium_current_Xs_gate (dimensionless).
 * RATES[7] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[8] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[9] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
 * RATES[10] is d/dt d in component L_type_Ca_current_d_gate (dimensionless).
 * RATES[11] is d/dt f in component L_type_Ca_current_f_gate (dimensionless).
 * RATES[12] is d/dt fCa in component L_type_Ca_current_fCa_gate (dimensionless).
 * RATES[13] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[14] is d/dt r in component transient_outward_current_r_gate (dimensionless).
 * RATES[16] is d/dt g in component calcium_dynamics (dimensionless).
 * RATES[3] is d/dt Ca_i in component calcium_dynamics (millimolar).
 * RATES[15] is d/dt Ca_SR in component calcium_dynamics (millimolar).
 * RATES[2] is d/dt Na_i in component sodium_dynamics (millimolar).
 * RATES[1] is d/dt K_i in component potassium_dynamics (millimolar).
 *
 * Adding 3 constant for multiple stimuli support:
 * CONSTANTS[46] is nstim1 in component membrane (dimensionless).
 * CONSTANTS[47] is stim_period_S2 in component membrane (millisecond).
 * CONSTANTS[48] is nstim2 in component membrane (dimensionless).
 * Adding a Algebraic variable for Iion output:
 * ALGEBRAIC[69] is Iion in component membrane (picoA_per_picoF).
 *
 */

#include "dictionary.H"
#include "scalar.H"
#include "TNNP_2004Names.H"


 static const char* TNNP_STATES_NAMES[17] = {
    "V",        // 0: membrane voltage (millivolt)
    "K_i",      // 1: intracellular potassium (millimolar)
    "Na_i",     // 2: intracellular sodium (millimolar)
    "Ca_i",     // 3: intracellular calcium (millimolar)
    "Xr1",      // 4: Xr1 gate (dimensionless)
    "Xr2",      // 5: Xr2 gate (dimensionless)
    "Xs",       // 6: Xs gate (dimensionless)
    "m",        // 7: m gate (dimensionless)
    "h",        // 8: h gate (dimensionless)
    "j",        // 9: j gate (dimensionless)
    "d",        // 10: d gate (dimensionless)
    "f",        // 11: f gate (dimensionless)
    "fCa",      // 12: fCa gate (dimensionless)
    "s",        // 13: s gate (dimensionless)
    "r",        // 14: r gate (dimensionless)
    "Ca_SR",    // 15: sarcoplasmic reticulum calcium (millimolar)
    "g"         // 16: g gate (dimensionless)
};
static const char* TNNP_ALGEBRAIC_NAMES[70] = {
    "I_stim",                    // 0
    "xr1_inf",                  // 1
    "xr2_inf",                  // 2
    "xs_inf",                   // 3
    "m_inf",                    // 4
    "h_inf",                    // 5
    "j_inf",                    // 6
    "d_inf",                    // 7
    "f_inf",                    // 8
    "alpha_fCa",               // 9
    "s_inf",                    // 10
    "r_inf",                    // 11
    "g_inf",                    // 12
    "E_Na",                     // 13
    "alpha_xr1",               // 14
    "alpha_xr2",               // 15
    "alpha_xs",                // 16
    "alpha_m",                 // 17
    "alpha_h",                 // 18
    "alpha_j",                 // 19
    "alpha_d",                 // 20
    "tau_f",                    // 21
    "beta_fCa",                // 22
    "tau_s",                    // 23
    "tau_r",                    // 24
    "d_g",                      // 25
    "E_K",                      // 26
    "beta_xr1",                // 27
    "beta_xr2",                // 28
    "beta_xs",                 // 29
    "beta_m",                  // 30
    "beta_h",                  // 31
    "beta_j",                  // 32
    "gama_fCa",                // 33
    "beta_d",                  // 34
    "E_Ks",                     // 35
    "tau_xr1",                 // 36
    "tau_xr2",                 // 37
    "tau_xs",                  // 38
    "tau_m",                   // 39
    "tau_h",                   // 40
    "tau_j",                   // 41
    "gamma_d",                 // 42
    "fCa_inf",                 // 43
    "E_Ca",                     // 44
    "tau_d",                   // 45
    "d_fCa",                   // 46
    "alpha_K1",                // 47
    "beta_K1",                 // 48
    "xK1_inf",                 // 49
    "i_K1",                     // 50
    "i_Kr",                     // 51
    "i_Ks",                     // 52
    "i_Na",                     // 53
    "i_b_Na",                   // 54
    "i_CaL",                    // 55
    "i_b_Ca",                   // 56
    "i_to",                     // 57
    "i_NaK",                    // 58
    "i_NaCa",                   // 59
    "i_p_Ca",                   // 60
    "i_p_K",                    // 61
    "i_rel",                    // 62
    "i_up",                     // 63
    "i_leak",                   // 64
    "ddt_Ca_i_total",           // 65
    "ddt_Ca_sr_total",          // 66
    "f_JCa_i_free",             // 67
    "f_JCa_sr_free",            // 68
    "Iion_cm"                   // 69
};

void
TNNPinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)

{
STATES[0] = -86.2;
CONSTANTS[0] = 8.314;
CONSTANTS[1] = 310;
CONSTANTS[2] = 96.485;
CONSTANTS[3] = 185;
CONSTANTS[4] = 16404;

// CONSTANTS[5] = 20;
// CONSTANTS[6] = 500;
// CONSTANTS[7] = 1;
// CONSTANTS[8] = -52;
//-----

CONSTANTS[9] = 0.03;
CONSTANTS[10] = 5.4;
CONSTANTS[11] = 140;
STATES[1] = 138.3;
STATES[2] = 11.6;
CONSTANTS[12] = 2;
STATES[3] = 0.0002;
CONSTANTS[13] = 5.405;
CONSTANTS[14] = 0.096;
STATES[4] = 0;
STATES[5] = 1;
//Conditional tissue conductance Gks
CONSTANTS[15] = (tissueFlag == 2)
          ? 0.062
          : 0.245;

STATES[6] = 0;
CONSTANTS[16] = 14.838;
STATES[7] = 0;
STATES[8] = 0.75;
STATES[9] = 0.75;
CONSTANTS[17] = 0.00029;
CONSTANTS[18] = 0.175;
STATES[10] = 0;
STATES[11] = 1;
STATES[12] = 1;
CONSTANTS[19] = 0.000592;
//Conditional tissue conductance Gto
CONSTANTS[20] = (tissueFlag == 3)
          ? 0.073
          : 0.294;

STATES[13] = 1;
STATES[14] = 0;
CONSTANTS[21] = 1.362;
CONSTANTS[22] = 1;
CONSTANTS[23] = 40;
CONSTANTS[24] = 1000;
CONSTANTS[25] = 0.1;
CONSTANTS[26] = 2.5;
CONSTANTS[27] = 0.35;
CONSTANTS[28] = 1.38;
CONSTANTS[29] = 87.5;
CONSTANTS[30] = 0.825;
CONSTANTS[31] = 0.0005;
CONSTANTS[32] = 0.0146;
STATES[15] = 0.2;
STATES[16] = 1;
CONSTANTS[33] = 2;
CONSTANTS[34] = 0.016464;
CONSTANTS[35] = 0.25;
CONSTANTS[36] = 0.008232;
CONSTANTS[37] = 0.00025;
CONSTANTS[38] = 8e-5;
CONSTANTS[39] = 0.000425;
CONSTANTS[40] = 0.15;
CONSTANTS[41] = 0.001;
CONSTANTS[42] = 10;
CONSTANTS[43] = 0.3;
CONSTANTS[44] = 1094;
CONSTANTS[45] = 2.00000;

}
void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver)
{
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[21] =  1125.00*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[0])/10.0000));
RATES[11] = (ALGEBRAIC[8] - STATES[11])/ALGEBRAIC[21];

//Conditional tissue flag for the Ito current
ALGEBRAIC[10] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + exp((STATES[0] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + exp((STATES[0] + 28.0000) / 6.00000));

ALGEBRAIC[23] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[0]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[0]+67.0000, 2.00000)/1000.000)+8.00000;

RATES[13] = (ALGEBRAIC[10] - STATES[13])/ALGEBRAIC[23];
ALGEBRAIC[11] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/6.00000));
ALGEBRAIC[24] =  9.50000*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+0.800000;
RATES[14] = (ALGEBRAIC[11] - STATES[14])/ALGEBRAIC[24];
ALGEBRAIC[12] = (STATES[3]<0.000350000 ? 1.00000/(1.00000+pow(STATES[3]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[3]/0.000350000, 16.0000)));
ALGEBRAIC[25] = (ALGEBRAIC[12] - STATES[16])/CONSTANTS[33];
RATES[16] = (ALGEBRAIC[12]>STATES[16]&&STATES[0]>- 60.0000 ? 0.00000 : ALGEBRAIC[25]);
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[14] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[27] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[36] =  1.00000*ALGEBRAIC[14]*ALGEBRAIC[27];
RATES[4] = (ALGEBRAIC[1] - STATES[4])/ALGEBRAIC[36];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[15] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[28] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[37] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[28];
RATES[5] = (ALGEBRAIC[2] - STATES[5])/ALGEBRAIC[37];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[16] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[29] = 1.00000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[29];
RATES[6] = (ALGEBRAIC[3] - STATES[6])/ALGEBRAIC[38];
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[17] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[30] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[30];
RATES[7] = (ALGEBRAIC[4] - STATES[7])/ALGEBRAIC[39];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[18] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[31] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[40] = 1.00000/(ALGEBRAIC[18]+ALGEBRAIC[31]);
RATES[8] = (ALGEBRAIC[5] - STATES[8])/ALGEBRAIC[40];
ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[32] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[41] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[32]);
RATES[9] = (ALGEBRAIC[6] - STATES[9])/ALGEBRAIC[41];
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/7.50000));
ALGEBRAIC[20] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[33] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[42] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[45] =  1.00000*ALGEBRAIC[20]*ALGEBRAIC[33]+ALGEBRAIC[42];
RATES[10] = (ALGEBRAIC[7] - STATES[10])/ALGEBRAIC[45];
ALGEBRAIC[9] = 1.00000/(1.00000+pow(STATES[3]/0.000325000, 8.00000));
ALGEBRAIC[22] = 0.100000/(1.00000+exp((STATES[3] - 0.000500000)/0.000100000));
ALGEBRAIC[34] = 0.200000/(1.00000+exp((STATES[3] - 0.000750000)/0.000800000));
ALGEBRAIC[43] = (ALGEBRAIC[9]+ALGEBRAIC[22]+ALGEBRAIC[34]+0.230000)/1.46000;
ALGEBRAIC[46] = (ALGEBRAIC[43] - STATES[12])/CONSTANTS[45];
RATES[12] = (ALGEBRAIC[43]>STATES[12]&&STATES[0]>- 60.0000 ? 0.00000 : ALGEBRAIC[46]);
ALGEBRAIC[58] = (( (( CONSTANTS[21]*CONSTANTS[10])/(CONSTANTS[10]+CONSTANTS[22]))*STATES[2])/(STATES[2]+CONSTANTS[23]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[11]/STATES[2]);
ALGEBRAIC[53] =  CONSTANTS[16]*pow(STATES[7], 3.00000)*STATES[8]*STATES[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[54] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[59] = ( CONSTANTS[24]*( exp(( CONSTANTS[27]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[12] -  exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[11], 3.00000)*STATES[3]*CONSTANTS[26]))/( (pow(CONSTANTS[29], 3.00000)+pow(CONSTANTS[11], 3.00000))*(CONSTANTS[28]+CONSTANTS[12])*(1.00000+ CONSTANTS[25]*exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
RATES[2] = ( - (ALGEBRAIC[53]+ALGEBRAIC[54]+ 3.00000*ALGEBRAIC[58]+ 3.00000*ALGEBRAIC[59])*CONSTANTS[3])/( CONSTANTS[4]*CONSTANTS[2]);
ALGEBRAIC[26] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/STATES[1]);
ALGEBRAIC[47] = 0.100000/(1.00000+exp( 0.0600000*((STATES[0] - ALGEBRAIC[26]) - 200.000)));
ALGEBRAIC[48] = ( 3.00000*exp( 0.000200000*((STATES[0] - ALGEBRAIC[26])+100.000))+ 1.00000*exp( 0.100000*((STATES[0] - ALGEBRAIC[26]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[0] - ALGEBRAIC[26])));
ALGEBRAIC[49] = ALGEBRAIC[47]/(ALGEBRAIC[47]+ALGEBRAIC[48]);
ALGEBRAIC[50] =  CONSTANTS[13]*ALGEBRAIC[49]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*(STATES[0] - ALGEBRAIC[26]);

//Conditional tissue conductance Ito
ALGEBRAIC[57] =  CONSTANTS[20]*STATES[14]*STATES[13]*(STATES[0] - ALGEBRAIC[26]);

ALGEBRAIC[51] =  CONSTANTS[14]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*STATES[4]*STATES[5]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[10]+ CONSTANTS[9]*CONSTANTS[11])/(STATES[1]+ CONSTANTS[9]*STATES[2]));
ALGEBRAIC[52] =  CONSTANTS[15]*pow(STATES[6], 2.00000)*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[55] = ( (( CONSTANTS[18]*STATES[10]*STATES[11]*STATES[12]*4.00000*STATES[0]*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( STATES[3]*exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) -  0.341000*CONSTANTS[12]))/(exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[44] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[12]/STATES[3]);
ALGEBRAIC[56] =  CONSTANTS[19]*(STATES[0] - ALGEBRAIC[44]);
ALGEBRAIC[61] = ( CONSTANTS[32]*(STATES[0] - ALGEBRAIC[26]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[60] = ( CONSTANTS[30]*STATES[3])/(STATES[3]+CONSTANTS[31]);

ALGEBRAIC[69] = ALGEBRAIC[50]+ALGEBRAIC[57]+ALGEBRAIC[51]+
                ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[58]+
                ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[59]+
                ALGEBRAIC[56]+ALGEBRAIC[61]+ALGEBRAIC[60];

if (solveVmWithinODESolver)
{

    double I_stim = 0.0;

    // --- S1 window ---
    double phase_S1 = VOI - floor(VOI/CONSTANTS[6])*CONSTANTS[6];
    bool in_S1 =
    (
        phase_S1 >= CONSTANTS[5] &&
        phase_S1 <= CONSTANTS[5] + CONSTANTS[7] &&
        VOI <= CONSTANTS[5] + CONSTANTS[6]*CONSTANTS[46]
    );
    if (in_S1)
    {
        I_stim = CONSTANTS[8];
    }
    else
    {
        // --- S2 window (ONLY if valid) ---
        double T2  = CONSTANTS[47];
        int    n2  = CONSTANTS[48];

        if (T2 > 0.0 && n2 > 0)
        {
            double tS1End = CONSTANTS[5]
                        + CONSTANTS[6]*CONSTANTS[46];

            bool in_S2_total_window =
            (
                VOI >= tS1End &&
                VOI <= tS1End + n2*T2
            );

            if (in_S2_total_window)
            {
                double phase = fmod(VOI - tS1End, T2);

                if (phase <= CONSTANTS[7])
                {
                    I_stim = CONSTANTS[8];
                }
            }
        }
    }

    ALGEBRAIC[0] = I_stim;

    RATES[0] = -(ALGEBRAIC[69]+ALGEBRAIC[0]);
}
else
{
    // Keep Vm constant: zero derivative
    RATES[0] = 0.0;
}

RATES[1] = ( - ((ALGEBRAIC[50]+ALGEBRAIC[57]+ALGEBRAIC[51]+ALGEBRAIC[52]+ALGEBRAIC[61]+ALGEBRAIC[0]) -  2.00000*ALGEBRAIC[58])*CONSTANTS[3])/( CONSTANTS[4]*CONSTANTS[2]);
ALGEBRAIC[62] =  (( CONSTANTS[34]*pow(STATES[15], 2.00000))/(pow(CONSTANTS[35], 2.00000)+pow(STATES[15], 2.00000))+CONSTANTS[36])*STATES[10]*STATES[16];
ALGEBRAIC[63] = CONSTANTS[39]/(1.00000+pow(CONSTANTS[37], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[64] =  CONSTANTS[38]*(STATES[15] - STATES[3]);
ALGEBRAIC[65] = (( (- ((ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[60]) -  2.00000*ALGEBRAIC[59])/( 2.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3]+ALGEBRAIC[64]) - ALGEBRAIC[63])+ALGEBRAIC[62];
ALGEBRAIC[67] = 1.00000/(1.00000+( CONSTANTS[40]*CONSTANTS[41])/pow(STATES[3]+CONSTANTS[41], 2.00000));
RATES[3] =  ALGEBRAIC[65]*ALGEBRAIC[67];
ALGEBRAIC[66] =  (CONSTANTS[4]/CONSTANTS[44])*(ALGEBRAIC[63] - (ALGEBRAIC[62]+ALGEBRAIC[64]));
ALGEBRAIC[68] = 1.00000/(1.00000+( CONSTANTS[42]*CONSTANTS[43])/pow(STATES[15]+CONSTANTS[43], 2.00000));
RATES[15] =  ALGEBRAIC[66]*ALGEBRAIC[68];
}

void
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{

ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+20.0000)/7.00000));
ALGEBRAIC[21] =  1125.00*exp(- pow(STATES[0]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[0])/10.0000));

//Conditional tissue flag for the Ito current
ALGEBRAIC[10] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + exp((STATES[0] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + exp((STATES[0] + 28.0000) / 6.00000));

ALGEBRAIC[23] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[0]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[0] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[0]+67.0000, 2.00000)/1000.000)+8.00000;

ALGEBRAIC[11] = 1.00000/(1.00000+exp((20.0000 - STATES[0])/6.00000));
ALGEBRAIC[24] =  9.50000*exp(- pow(STATES[0]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[12] = (STATES[3]<0.000350000 ? 1.00000/(1.00000+pow(STATES[3]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[3]/0.000350000, 16.0000)));
ALGEBRAIC[25] = (ALGEBRAIC[12] - STATES[16])/CONSTANTS[33];
ALGEBRAIC[1] = 1.00000/(1.00000+exp((- 26.0000 - STATES[0])/7.00000));
ALGEBRAIC[14] = 450.000/(1.00000+exp((- 45.0000 - STATES[0])/10.0000));
ALGEBRAIC[27] = 6.00000/(1.00000+exp((STATES[0]+30.0000)/11.5000));
ALGEBRAIC[36] =  1.00000*ALGEBRAIC[14]*ALGEBRAIC[27];
ALGEBRAIC[2] = 1.00000/(1.00000+exp((STATES[0]+88.0000)/24.0000));
ALGEBRAIC[15] = 3.00000/(1.00000+exp((- 60.0000 - STATES[0])/20.0000));
ALGEBRAIC[28] = 1.12000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[37] =  1.00000*ALGEBRAIC[15]*ALGEBRAIC[28];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/14.0000));
ALGEBRAIC[16] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[0])/6.00000)), 1.0 / 2);
ALGEBRAIC[29] = 1.00000/(1.00000+exp((STATES[0] - 60.0000)/20.0000));
ALGEBRAIC[38] =  1.00000*ALGEBRAIC[16]*ALGEBRAIC[29];
ALGEBRAIC[4] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[0])/9.03000), 2.00000);
ALGEBRAIC[17] = 1.00000/(1.00000+exp((- 60.0000 - STATES[0])/5.00000));
ALGEBRAIC[30] = 0.100000/(1.00000+exp((STATES[0]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[0] - 50.0000)/200.000));
ALGEBRAIC[39] =  1.00000*ALGEBRAIC[17]*ALGEBRAIC[30];
ALGEBRAIC[5] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[18] = (STATES[0]<- 40.0000 ?  0.0570000*exp(- (STATES[0]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[31] = (STATES[0]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[0])+ 310000.*exp( 0.348500*STATES[0]) : 0.770000/( 0.130000*(1.00000+exp((STATES[0]+10.6600)/- 11.1000))));
ALGEBRAIC[40] = 1.00000/(ALGEBRAIC[18]+ALGEBRAIC[31]);
ALGEBRAIC[6] = 1.00000/pow(1.00000+exp((STATES[0]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[19] = (STATES[0]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[0]) -  6.94800e-06*exp( - 0.0439100*STATES[0]))*(STATES[0]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[0]+79.2300))) : 0.00000);
ALGEBRAIC[32] = (STATES[0]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[0]))/(1.00000+exp( - 0.137800*(STATES[0]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[0]))/(1.00000+exp( - 0.100000*(STATES[0]+32.0000))));
ALGEBRAIC[41] = 1.00000/(ALGEBRAIC[19]+ALGEBRAIC[32]);
ALGEBRAIC[7] = 1.00000/(1.00000+exp((- 5.00000 - STATES[0])/7.50000));
ALGEBRAIC[20] = 1.40000/(1.00000+exp((- 35.0000 - STATES[0])/13.0000))+0.250000;
ALGEBRAIC[33] = 1.40000/(1.00000+exp((STATES[0]+5.00000)/5.00000));
ALGEBRAIC[42] = 1.00000/(1.00000+exp((50.0000 - STATES[0])/20.0000));
ALGEBRAIC[45] =  1.00000*ALGEBRAIC[20]*ALGEBRAIC[33]+ALGEBRAIC[42];
ALGEBRAIC[9] = 1.00000/(1.00000+pow(STATES[3]/0.000325000, 8.00000));
ALGEBRAIC[22] = 0.100000/(1.00000+exp((STATES[3] - 0.000500000)/0.000100000));
ALGEBRAIC[34] = 0.200000/(1.00000+exp((STATES[3] - 0.000750000)/0.000800000));
ALGEBRAIC[43] = (ALGEBRAIC[9]+ALGEBRAIC[22]+ALGEBRAIC[34]+0.230000)/1.46000;
ALGEBRAIC[46] = (ALGEBRAIC[43] - STATES[12])/CONSTANTS[45];
ALGEBRAIC[58] = (( (( CONSTANTS[21]*CONSTANTS[10])/(CONSTANTS[10]+CONSTANTS[22]))*STATES[2])/(STATES[2]+CONSTANTS[23]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))+ 0.0353000*exp(( - STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[13] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[11]/STATES[2]);
ALGEBRAIC[53] =  CONSTANTS[16]*pow(STATES[7], 3.00000)*STATES[8]*STATES[9]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[54] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[13]);
ALGEBRAIC[59] = ( CONSTANTS[24]*( exp(( CONSTANTS[27]*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(STATES[2], 3.00000)*CONSTANTS[12] -  exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*pow(CONSTANTS[11], 3.00000)*STATES[3]*CONSTANTS[26]))/( (pow(CONSTANTS[29], 3.00000)+pow(CONSTANTS[11], 3.00000))*(CONSTANTS[28]+CONSTANTS[12])*(1.00000+ CONSTANTS[25]*exp(( (CONSTANTS[27] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))));
ALGEBRAIC[26] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[10]/STATES[1]);
ALGEBRAIC[47] = 0.100000/(1.00000+exp( 0.0600000*((STATES[0] - ALGEBRAIC[26]) - 200.000)));
ALGEBRAIC[48] = ( 3.00000*exp( 0.000200000*((STATES[0] - ALGEBRAIC[26])+100.000))+ 1.00000*exp( 0.100000*((STATES[0] - ALGEBRAIC[26]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[0] - ALGEBRAIC[26])));
ALGEBRAIC[49] = ALGEBRAIC[47]/(ALGEBRAIC[47]+ALGEBRAIC[48]);
ALGEBRAIC[50] =  CONSTANTS[13]*ALGEBRAIC[49]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[57] =  CONSTANTS[20]*STATES[14]*STATES[13]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[51] =  CONSTANTS[14]* pow((CONSTANTS[10]/5.40000), 1.0 / 2)*STATES[4]*STATES[5]*(STATES[0] - ALGEBRAIC[26]);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log((CONSTANTS[10]+ CONSTANTS[9]*CONSTANTS[11])/(STATES[1]+ CONSTANTS[9]*STATES[2]));
ALGEBRAIC[52] =  CONSTANTS[15]*pow(STATES[6], 2.00000)*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[55] = ( (( CONSTANTS[18]*STATES[10]*STATES[11]*STATES[12]*4.00000*STATES[0]*pow(CONSTANTS[2], 2.00000))/( CONSTANTS[0]*CONSTANTS[1]))*( STATES[3]*exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) -  0.341000*CONSTANTS[12]))/(exp(( 2.00000*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[44] =  (( 0.500000*CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(CONSTANTS[12]/STATES[3]);
ALGEBRAIC[56] =  CONSTANTS[19]*(STATES[0] - ALGEBRAIC[44]);
ALGEBRAIC[61] = ( CONSTANTS[32]*(STATES[0] - ALGEBRAIC[26]))/(1.00000+exp((25.0000 - STATES[0])/5.98000));
ALGEBRAIC[60] = ( CONSTANTS[30]*STATES[3])/(STATES[3]+CONSTANTS[31]);

ALGEBRAIC[62] =  (( CONSTANTS[34]*pow(STATES[15], 2.00000))/(pow(CONSTANTS[35], 2.00000)+pow(STATES[15], 2.00000))+CONSTANTS[36])*STATES[10]*STATES[16];
ALGEBRAIC[63] = CONSTANTS[39]/(1.00000+pow(CONSTANTS[37], 2.00000)/pow(STATES[3], 2.00000));
ALGEBRAIC[64] =  CONSTANTS[38]*(STATES[15] - STATES[3]);
ALGEBRAIC[65] = (( (- ((ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[60]) -  2.00000*ALGEBRAIC[59])/( 2.00000*CONSTANTS[4]*CONSTANTS[2]))*CONSTANTS[3]+ALGEBRAIC[64]) - ALGEBRAIC[63])+ALGEBRAIC[62];
ALGEBRAIC[67] = 1.00000/(1.00000+( CONSTANTS[40]*CONSTANTS[41])/pow(STATES[3]+CONSTANTS[41], 2.00000));
ALGEBRAIC[66] =  (CONSTANTS[4]/CONSTANTS[44])*(ALGEBRAIC[63] - (ALGEBRAIC[62]+ALGEBRAIC[64]));
ALGEBRAIC[68] = 1.00000/(1.00000+( CONSTANTS[42]*CONSTANTS[43])/pow(STATES[15]+CONSTANTS[43], 2.00000));

ALGEBRAIC[69] = ALGEBRAIC[50]+ALGEBRAIC[57]+ALGEBRAIC[51]+
                ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[58]+
                ALGEBRAIC[53]+ALGEBRAIC[54]+ALGEBRAIC[59]+
                ALGEBRAIC[56]+ALGEBRAIC[61]+ALGEBRAIC[60];


}
