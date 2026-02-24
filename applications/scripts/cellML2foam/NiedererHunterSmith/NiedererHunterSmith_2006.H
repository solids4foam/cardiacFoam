#include "NiedererHunterSmith_2006Names.H"

static const char* NiedererHunterSmithSTATES_NAMES[NUM_STATES] = {
    "Ca_i_TRPN",
    "tropomyosin_z",
    "Cross_Bridges_Q1",
    "Cross_Bridges_Q2",
    "Cross_Bridges_Q3",
};

static const char* NiedererHunterSmithALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_ExtensionRatio",
    "AV_lambda",
    "AV_lambda_prev",
    "AV_overlap",
    "AV_Q",
    "AV_Ca_i",
    "AV_Ca_50",
    "AV_beta_Tm",
    "AV_Ca_b",
    "AV_Ca_TRPN_50",
    "AV_alpha_Tm",
    "AV_z_max",
    "AV_T_Base",
    "AV_T_0",
    "AV_Tension",
    "AV_k_off",
    "AV_J_TRPN",
};


static inline const Foam::HashTable<Foam::wordList>&
NiedererHunterSmithDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}
void
NiedererHunterSmithinitConsts(double* CONSTANTS, double* RATES, double* STATES)
{
/* Myofilaments */
    CONSTANTS[AC_dExtensionRatiodt] = 0.0;
    
    /* filament_overlap */
    CONSTANTS[AC_beta_0] = 4.9;
    
    /* Cross_Bridges */
    CONSTANTS[AC_A_1] = -29.0;
    CONSTANTS[AC_A_2] = 138.0;
    CONSTANTS[AC_A_3] = 129.0;
    CONSTANTS[AC_a] = 0.35;
    CONSTANTS[AC_alpha_1] = 0.03;
    CONSTANTS[AC_alpha_2] = 0.13;
    CONSTANTS[AC_alpha_3] = 0.625;
    
    /* length_independent_tension */
    CONSTANTS[AC_T_ref] = 56.2;
    
    /* tropomyosin */
    CONSTANTS[AC_Ca_50ref] = 0.00105;
    CONSTANTS[AC_K_z] = 0.15;
    CONSTANTS[AC_alpha_0] = 0.008;
    CONSTANTS[AC_alpha_r1] = 0.002;
    CONSTANTS[AC_alpha_r2] = 0.00175;
    CONSTANTS[AC_beta_1] = -4.0;
    CONSTANTS[AC_n_Hill] = 3.0;
    CONSTANTS[AC_n_Rel] = 3.0;
    CONSTANTS[AC_z_p] = 0.85;
    CONSTANTS[AC_K_1] = CONSTANTS[AC_alpha_r2] * pow(CONSTANTS[AC_z_p],
                                                     CONSTANTS[AC_n_Rel] - 1.0) * CONSTANTS[AC_n_Rel] * pow(CONSTANTS[AC_K_z],
                                                                                                            CONSTANTS[AC_n_Rel]) / pow(pow(CONSTANTS[AC_z_p], CONSTANTS[AC_n_Rel]) + pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel]),
                                                                                                                                       2.0);
    CONSTANTS[AC_K_2] = CONSTANTS[AC_alpha_r2] * pow(CONSTANTS[AC_z_p],
                                                     CONSTANTS[AC_n_Rel]) / (pow(CONSTANTS[AC_z_p], CONSTANTS[AC_n_Rel]) + pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel])) * (1.0 - CONSTANTS[AC_n_Rel] * pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel]) / (pow(CONSTANTS[AC_z_p], CONSTANTS[AC_n_Rel]) + pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel])));
    
    /* troponin */
    CONSTANTS[AC_Ca_TRPN_Max] = 0.07;
    CONSTANTS[AC_gamma_trpn] = 2.0;
    CONSTANTS[AC_k_Ref_off] = 0.2;
    CONSTANTS[AC_k_on] = 100.0;

    /* --- Initial values --- */

    STATES[Ca_i_TRPN] =  6.75931398649999987e-02;
    STATES[tropomyosin_z] =  1.44179378369999993e-02;
    STATES[Cross_Bridges_Q1] = 0.0;
    STATES[Cross_Bridges_Q2] = 0.0;
    STATES[Cross_Bridges_Q3] = 0.0;
}


void
NiedererHunterSmithcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver)
{
    
    /* Myofilaments */
    ALGEBRAIC[AV_ExtensionRatio] = ((VOI > 300000.0) ? 1.0 : 1.0);
    ALGEBRAIC[AV_lambda] = (((ALGEBRAIC[AV_ExtensionRatio] > 0.8) && (ALGEBRAIC[AV_ExtensionRatio] <= 1.15)) ? ALGEBRAIC[AV_ExtensionRatio] : ((ALGEBRAIC[AV_ExtensionRatio] > 1.15) ? 1.15 : 0.8));
    ALGEBRAIC[AV_lambda_prev] = ALGEBRAIC[AV_ExtensionRatio];
    
    /* filament_overlap */
    ALGEBRAIC[AV_overlap] = 1.0 + CONSTANTS[AC_beta_0] * (ALGEBRAIC[AV_lambda] - 1.0);
    
    /* Cross_Bridges */
    ALGEBRAIC[AV_Q] = STATES[Cross_Bridges_Q1] + STATES[Cross_Bridges_Q2] + STATES[Cross_Bridges_Q3];
    RATES[Cross_Bridges_Q1] = CONSTANTS[AC_A_1] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_1] * STATES[Cross_Bridges_Q1];
    RATES[Cross_Bridges_Q2] = CONSTANTS[AC_A_2] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_2] * STATES[Cross_Bridges_Q2];
    RATES[Cross_Bridges_Q3] = CONSTANTS[AC_A_3] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_3] * STATES[Cross_Bridges_Q3];
    
    /* intracellular_ion_concentrations */
    ALGEBRAIC[AV_Ca_i] = ((VOI < 1.0) ? 1000.0 * 1.84330000000000010e-07 : (((VOI >= 10.0) && (VOI < 15.0)) ? 1000.0 * (1.055 * pow(VOI / 1000.0, 3.0) - 0.03507 * pow(VOI / 1000.0, 2.0) + 0.0003992 * VOI / 1000.0 - 1.356e-06) : (((VOI >= 15.0) && (VOI < 55.0)) ? 1000.0 * (0.014 * pow(VOI / 1000.0, 3.0) - 0.002555 * pow(VOI / 1000.0, 2.0) + 0.0001494 * VOI / 1000.0 - 1.428e-06) : (((VOI >= 55.0) && (VOI < 250.0)) ? 1000.0 * (1.739e-05 * pow(VOI / 1000.0, 3.0) - 3.209e-06 * pow(VOI / 1000.0, 2.0) - 5.689e-06 * VOI / 1000.0 + 1.719e-06) : (((VOI >= 250.0) && (VOI < 490.0)) ? 1000.0 * (0.0001321 * pow(VOI / 1000.0, 4.0) - 0.0002197 * pow(VOI / 1000.0, 3.0) + 0.0001374 * pow(VOI / 1000.0, 2.0) - 3.895e-05 * VOI / 1000.0 + 4.441e-06) : 1000.0 * 1.21480000000000000e-07)))));
    
    /* tropomyosin */
    ALGEBRAIC[AV_Ca_50] = CONSTANTS[AC_Ca_50ref] * (1.0 + CONSTANTS[AC_beta_1] * (ALGEBRAIC[AV_lambda] - 1.0));
    ALGEBRAIC[AV_beta_Tm] = CONSTANTS[AC_alpha_r1] + CONSTANTS[AC_alpha_r2] * pow(STATES[tropomyosin_z],
                                                                       CONSTANTS[AC_n_Rel] - 1.0) / (pow(STATES[tropomyosin_z], CONSTANTS[AC_n_Rel]) + pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel]));
    
    /* *remaining* */
    ALGEBRAIC[AV_Ca_b] = CONSTANTS[AC_Ca_TRPN_Max] - STATES[Ca_i_TRPN];
    ALGEBRAIC[AV_Ca_TRPN_50] = ALGEBRAIC[AV_Ca_50] * CONSTANTS[AC_Ca_TRPN_Max] / (ALGEBRAIC[AV_Ca_50] + CONSTANTS[AC_k_Ref_off] / CONSTANTS[AC_k_on] * (1.0 - (1.0 + CONSTANTS[AC_beta_0] * (ALGEBRAIC[AV_lambda] - 1.0)) * 0.5 / CONSTANTS[AC_gamma_trpn]));
    ALGEBRAIC[AV_alpha_Tm] = CONSTANTS[AC_alpha_0] * pow(ALGEBRAIC[AV_Ca_b] / ALGEBRAIC[AV_Ca_TRPN_50],
                                                         CONSTANTS[AC_n_Hill]);
    ALGEBRAIC[AV_z_max] = (CONSTANTS[AC_alpha_0] / pow(ALGEBRAIC[AV_Ca_TRPN_50] / CONSTANTS[AC_Ca_TRPN_Max], CONSTANTS[AC_n_Hill]) - CONSTANTS[AC_K_2]) / (CONSTANTS[AC_alpha_r1] + CONSTANTS[AC_K_1] + CONSTANTS[AC_alpha_0] / pow(ALGEBRAIC[AV_Ca_TRPN_50] / CONSTANTS[AC_Ca_TRPN_Max], CONSTANTS[AC_n_Hill]));
    ALGEBRAIC[AV_T_Base] = CONSTANTS[AC_T_ref] * STATES[tropomyosin_z] / ALGEBRAIC[AV_z_max];
    RATES[tropomyosin_z] = ALGEBRAIC[AV_alpha_Tm] * (1.0 - STATES[tropomyosin_z]) - ALGEBRAIC[AV_beta_Tm] * STATES[tropomyosin_z];
    ALGEBRAIC[AV_T_0] = ALGEBRAIC[AV_T_Base] * ALGEBRAIC[AV_overlap];
    ALGEBRAIC[AV_Tension] = ((ALGEBRAIC[AV_Q] < 0.0) ? ALGEBRAIC[AV_T_0] * (CONSTANTS[AC_a] * ALGEBRAIC[AV_Q] + 1.0) / (1.0 - ALGEBRAIC[AV_Q]) : ALGEBRAIC[AV_T_0] * (1.0 + (CONSTANTS[AC_a] + 2.0) * ALGEBRAIC[AV_Q]) / (1.0 + ALGEBRAIC[AV_Q]));
    ALGEBRAIC[AV_k_off] = ((1.0 - ALGEBRAIC[AV_Tension] / (CONSTANTS[AC_gamma_trpn] * CONSTANTS[AC_T_ref]) > 0.1) ? CONSTANTS[AC_k_Ref_off] * (1.0 - ALGEBRAIC[AV_Tension] / (CONSTANTS[AC_gamma_trpn] * CONSTANTS[AC_T_ref])) : CONSTANTS[AC_k_Ref_off] * 0.1);
    ALGEBRAIC[AV_J_TRPN] = (CONSTANTS[AC_Ca_TRPN_Max] - STATES[Ca_i_TRPN]) * ALGEBRAIC[AV_k_off] - ALGEBRAIC[AV_Ca_i] * STATES[Ca_i_TRPN] * CONSTANTS[AC_k_on];
    RATES[Ca_i_TRPN] = ALGEBRAIC[AV_J_TRPN];
    

    
}
