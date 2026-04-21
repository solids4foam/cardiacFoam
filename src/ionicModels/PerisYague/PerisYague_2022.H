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

#include "PerisYague_2022Names.H"
#include "stimulusIO.H"
#include "dictionary.H"
#include <cmath>

const char* PerisYague_2022STATES_NAMES[NUM_STATES] = {
    "membrane_V",
    "sodium_Nai",
    "potassium_Ki",
    "chloride_Cli",
    "calcium_Cai",
    "calcium_CaUp",
    "calcium_CaRel",
    "ina_m",
    "ina_h",
    "ina_j",
    "ikur_ua",
    "ikur_uif",
    "ikur_uis",
    "ikr_xr",
    "iks_xs",
    "ical_d",
    "ical_f",
    "ical_fCa",
    "iclca_qCa",
    "ryr_u",
    "ryr_v",
    "ryr_w"
};

const char* PerisYague_2022ALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_ICaL",
    "AV_IK1",
    "AV_IKr",
    "AV_IKs",
    "AV_IKur",
    "AV_INa",
    "AV_INaCa",
    "AV_INaK",
    "AV_IClCa",
    "AV_IbCa",
    "AV_IbNa",
    "AV_IbK",
    "AV_IpCa",

    "AV_ina_m_inf",
    "AV_ina_m_tau",
    "AV_ina_h_inf",
    "AV_ina_h_tau",
    "AV_ina_j_inf",
    "AV_ina_j_tau",

    "AV_ical_d_inf",
    "AV_ical_d_tau",
    "AV_ical_f_inf",
    "AV_ical_f_tau",

    "AV_ikur_ua_inf",
    "AV_ikur_ua_tau",
    "AV_ikur_uif_inf",
    "AV_ikur_uif_tau",
    "AV_ikur_uis_inf",
    "AV_ikur_uis_tau",
    "AV_gKur",

    "AV_ikr_xr_inf",
    "AV_ikr_xr_tau",

    "AV_iks_xs_inf",
    "AV_iks_xs_tau",

    "AV_iclca_qCa_inf",
    "AV_iclca_qCa_tau",

    "AV_fNaK",

    "AV_ryr_w_inf",
    "AV_ryr_w_tau",

    "Iion_cm",
    "Istim"
};

const char* PerisYague_2022CONSTANTS_NAMES[NUM_CONSTANTS] = {
    "AC_CMDN_max",
    "AC_CSQN_max",
    "AC_Ca_up_max",
    "AC_Cao",
    "AC_Cm",
    "AC_ECaL",
    "AC_F",
    "AC_FRT",
    "AC_INaCa_max",
    "AC_INaK_max",
    "AC_I_diff",
    "AC_I_up_max",
    "AC_KQ10",
    "AC_krel",
    "AC_K_up",
    "AC_KmCa",
    "AC_KmKo",
    "AC_KmNa",
    "AC_KmNai",
    "AC_Km_CMDN",
    "AC_Km_CSQN",
    "AC_Km_TRPN",
    "AC_Ko",
    "AC_Nao",
    "AC_Clo",
    "AC_R",
    "AC_RTF",
    "AC_T",
    "AC_TRPN_max",
    "AC_V_cell",
    "AC_V_i",
    "AC_V_rel",
    "AC_V_up",
    "AC_c1",
    "AC_c2",
    "AC_cajsr_u_tau",
    "AC_gNa",
    "AC_gK1",
    "AC_gKr",
    "AC_gKs",
    "AC_gKur_amp",
    "AC_gClCa",
    "AC_gbCa",
    "AC_gbNa",
    "AC_gbK",
    "AC_ical_fCa_tau",
    "AC_gamma",
    "AC_ksat",
    "AC_sigma",
    "AC_tau_tr",
    "AC_gCaL",
    "AC_IpCa_max"
};

static inline const Foam::HashTable<Foam::wordList>&
PerisYague_2022DependencyMap()
{
    static const Foam::HashTable<Foam::wordList> dep
    (
        {
            { "INa", Foam::wordList({ "AV_INa", "ina_m", "ina_h", "ina_j", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau" }) },
            { "I_Na", Foam::wordList({ "AV_INa", "ina_m", "ina_h", "ina_j", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau" }) },

            { "IKr", Foam::wordList({ "AV_IKr", "ikr_xr", "AV_ikr_xr_inf", "AV_ikr_xr_tau" }) },
            { "I_Kr", Foam::wordList({ "AV_IKr", "ikr_xr", "AV_ikr_xr_inf", "AV_ikr_xr_tau" }) },

            { "IKs", Foam::wordList({ "AV_IKs", "iks_xs", "AV_iks_xs_inf", "AV_iks_xs_tau" }) },
            { "I_Ks", Foam::wordList({ "AV_IKs", "iks_xs", "AV_iks_xs_inf", "AV_iks_xs_tau" }) },

            { "IKur", Foam::wordList({ "AV_IKur", "ikur_ua", "ikur_uif", "ikur_uis", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_uif_inf", "AV_ikur_uif_tau", "AV_ikur_uis_inf", "AV_ikur_uis_tau", "AV_gKur" }) },
            { "I_Kur", Foam::wordList({ "AV_IKur", "ikur_ua", "ikur_uif", "ikur_uis", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_uif_inf", "AV_ikur_uif_tau", "AV_ikur_uis_inf", "AV_ikur_uis_tau", "AV_gKur" }) },

            { "ICaL", Foam::wordList({ "AV_ICaL", "ical_d", "ical_f", "ical_fCa", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau" }) },
            { "I_CaL", Foam::wordList({ "AV_ICaL", "ical_d", "ical_f", "ical_fCa", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau" }) },

            { "IClCa", Foam::wordList({ "AV_IClCa", "iclca_qCa", "AV_iclca_qCa_inf", "AV_iclca_qCa_tau" }) },
            { "I_ClCa", Foam::wordList({ "AV_IClCa", "iclca_qCa", "AV_iclca_qCa_inf", "AV_iclca_qCa_tau" }) },

            { "IK1", Foam::wordList({ "AV_IK1" }) },
            { "I_K1", Foam::wordList({ "AV_IK1" }) },

            { "INaK", Foam::wordList({ "AV_INaK", "AV_fNaK" }) },
            { "I_NaK", Foam::wordList({ "AV_INaK", "AV_fNaK" }) },

            { "INaCa", Foam::wordList({ "AV_INaCa" }) },
            { "I_NaCa", Foam::wordList({ "AV_INaCa" }) },

            { "Ib", Foam::wordList({ "AV_IbNa", "AV_IbCa", "AV_IbK" }) },
            { "I_b", Foam::wordList({ "AV_IbNa", "AV_IbCa", "AV_IbK" }) },

            { "Iion", Foam::wordList({ "AV_INa", "AV_IKr", "AV_IKs", "AV_IKur", "AV_ICaL", "AV_IClCa", "AV_IK1", "AV_INaK", "AV_INaCa", "AV_IbNa", "AV_IbCa", "AV_IbK", "AV_IpCa", "Iion_cm" }) },
            { "I_ion", Foam::wordList({ "AV_INa", "AV_IKr", "AV_IKs", "AV_IKur", "AV_ICaL", "AV_IClCa", "AV_IK1", "AV_INaK", "AV_INaCa", "AV_IbNa", "AV_IbCa", "AV_IbK", "AV_IpCa", "Iion_cm" }) }
        }
    );
    return dep;
}

inline Foam::scalar PerisYague_2022computeIstim
(
    Foam::scalar t,
    const Foam::StimulusProtocol& stimulus
)
{
    return Foam::stimulusIO::computeStimulus(t, stimulus);
}

void
PerisYague_2022initConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
    // Table S1
    CONSTANTS[AC_R] = 8.3143;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_F] = 96.4867;
    CONSTANTS[AC_Cm] = 100.0;
    CONSTANTS[AC_V_cell] = 20100.0;
    CONSTANTS[AC_V_i] = 13668.0;
    CONSTANTS[AC_V_up] = 1109.52;
    CONSTANTS[AC_V_rel] = 96.48;
    CONSTANTS[AC_Ko] = 5.4;
    CONSTANTS[AC_Nao] = 140.0;
    CONSTANTS[AC_Cao] = 1.8;
    CONSTANTS[AC_Clo] = 132.0;

    CONSTANTS[AC_gNa] = 13.99;
    CONSTANTS[AC_gK1] = 0.08218;
    CONSTANTS[AC_gClCa] = 0.15731;
    CONSTANTS[AC_gKur_amp] = 0.45539;
    CONSTANTS[AC_gKr] = 0.01730;
    CONSTANTS[AC_gKs] = 0.0594;
    CONSTANTS[AC_gCaL] = 0.06574;
    CONSTANTS[AC_gbCa] = 0.00113;
    CONSTANTS[AC_gbNa] = 0.000674;
    CONSTANTS[AC_gbK] = 0.0; // Minimal leak if not specified

    CONSTANTS[AC_INaK_max] = 0.94935;
    CONSTANTS[AC_INaCa_max] = 2304.0;
    CONSTANTS[AC_IpCa_max] = 0.275;
    CONSTANTS[AC_I_up_max] = 0.005;
    CONSTANTS[AC_KQ10] = 3.0;
    CONSTANTS[AC_gamma] = 0.35;
    CONSTANTS[AC_KmNai] = 10.0;
    CONSTANTS[AC_KmKo] = 1.5;
    CONSTANTS[AC_KmNa] = 87.5;
    CONSTANTS[AC_KmCa] = 1.38;
    CONSTANTS[AC_ksat] = 0.1;

    CONSTANTS[AC_krel] = 30.0;
    CONSTANTS[AC_K_up] = 0.00092;
    CONSTANTS[AC_Ca_up_max] = 15.0;
    CONSTANTS[AC_CMDN_max] = 0.05;
    CONSTANTS[AC_TRPN_max] = 0.07;
    CONSTANTS[AC_CSQN_max] = 10.0;
    CONSTANTS[AC_Km_CMDN] = 0.00238;
    CONSTANTS[AC_Km_TRPN] = 0.0005;
    CONSTANTS[AC_Km_CSQN] = 0.8;
    CONSTANTS[AC_tau_tr] = 180.0;

    CONSTANTS[AC_c1] = 1.367e-13; // S65: Fn - 1.367e-13
    CONSTANTS[AC_c2] = 13.67e-16; // S65: / 13.67e-16
    CONSTANTS[AC_cajsr_u_tau] = 8.0;
    CONSTANTS[AC_ical_fCa_tau] = 2.0;

    CONSTANTS[AC_RTF] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F];
    CONSTANTS[AC_FRT] = 1.0 / CONSTANTS[AC_RTF];
    CONSTANTS[AC_sigma] = (exp(CONSTANTS[AC_Nao] / 67.3) - 1.0) / 7.0;
    CONSTANTS[AC_ECaL] = 65.0;

    // Table S2 Initial Values
    STATES[membrane_V] = -76.0;
    STATES[sodium_Nai] = 11.17;
    STATES[potassium_Ki] = 139.0;
    STATES[chloride_Cli] = 30.0;
    STATES[calcium_Cai] = 0.0001013;
    STATES[calcium_CaUp] = 1.488;
    STATES[calcium_CaRel] = 1.488;
    STATES[ina_m] = 2.91e-3;
    STATES[ina_h] = 0.965;
    STATES[ina_j] = 0.978;
    STATES[ikur_ua] = 4.96e-3;
    STATES[ikur_uif] = 1.0;
    STATES[ikur_uis] = 1.0;
    STATES[ikr_xr] = 3.29e-5;
    STATES[iks_xs] = 1.87e-2;
    STATES[ical_d] = 1.37e-4;
    STATES[ical_f] = 0.999;
    STATES[ical_fCa] = 0.775;
    STATES[iclca_qCa] = 0.0;
    STATES[ryr_u] = 0.0235;
    STATES[ryr_v] = 1.0;
    STATES[ryr_w] = 0.999;
}

void
PerisYague_2022computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
{
    double V = STATES[membrane_V];
    double Cai = STATES[calcium_Cai];
    double Nai = STATES[sodium_Nai];
    double Ki = STATES[potassium_Ki];
    double Cli = STATES[chloride_Cli];
    double CaRel = STATES[calcium_CaRel];
    double CaUp = STATES[calcium_CaUp];

    // Nernst Potentials
    double ENa = CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Nao] / Nai);
    double EK = CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Ko] / Ki);
    double ECa = 0.5 * CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Cao] / Cai);
    double ECl = -CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Clo] / Cli);

    // INa
    double am = (fabs(V + 47.13) < 1e-6) ? 3.2 : 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
    double bm = 0.08 * exp(-V / 11.0);
    ALGEBRAIC[AV_ina_m_inf] = am / (am + bm);
    ALGEBRAIC[AV_ina_m_tau] = 1.7 / (am + bm);
    RATES[ina_m] = (ALGEBRAIC[AV_ina_m_inf] - STATES[ina_m]) / ALGEBRAIC[AV_ina_m_tau];

    double ah = (V < -40.0) ? 0.135 * exp(-(V + 80.0) / 6.8) : 0.0;
    double bh = (V < -40.0) ? 3.56 * exp(0.079 * V) + 3.1e5 * exp(0.35 * V) : 0.13 / (1.0 + exp(-(V + 10.66) / 11.1));
    ALGEBRAIC[AV_ina_h_inf] = ah / (ah + bh);
    ALGEBRAIC[AV_ina_h_tau] = 2.0 / (ah + bh);
    RATES[ina_h] = (ALGEBRAIC[AV_ina_h_inf] - STATES[ina_h]) / ALGEBRAIC[AV_ina_h_tau];

    double aj = (V < -40.0) ? (-127140.0 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) * (V + 37.78) / (1.0 + exp(0.311 * (V + 79.23))) : 0.0;
    double bj = (V < -40.0) ? 0.1212 * exp(-0.01052 * V) / (1.0 + exp(-0.1378 * (V + 40.14))) : 0.3 * exp(-2.535e-7 * V) / (1.0 + exp(-0.1 * (V + 32.0)));
    ALGEBRAIC[AV_ina_j_inf] = aj / (aj + bj);
    ALGEBRAIC[AV_ina_j_tau] = 2.0 / (aj + bj);
    RATES[ina_j] = (ALGEBRAIC[AV_ina_j_inf] - STATES[ina_j]) / ALGEBRAIC[AV_ina_j_tau];
    ALGEBRAIC[AV_INa] = CONSTANTS[AC_gNa] * pow(STATES[ina_m], 3.0) * STATES[ina_h] * STATES[ina_j] * (V - ENa);

    // IK1
    ALGEBRAIC[AV_IK1] = CONSTANTS[AC_gK1] * (V - EK - 5.0) / (1.0 + exp(0.063 * (V + 70.0)));

    // IKr
    double axr = 0.0003 * (V + 14.1) / (1.0 - exp(-(V + 14.1) / 5.0));
    double bxr = 7.3898e-5 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1.0);
    ALGEBRAIC[AV_ikr_xr_inf] = axr / (axr + bxr);
    ALGEBRAIC[AV_ikr_xr_tau] = 1.0 / (axr + bxr);
    RATES[ikr_xr] = (ALGEBRAIC[AV_ikr_xr_inf] - STATES[ikr_xr]) / ALGEBRAIC[AV_ikr_xr_tau];
    ALGEBRAIC[AV_IKr] = CONSTANTS[AC_gKr] * STATES[ikr_xr] * (V - EK) / (1.0 + exp((V + 15.0) / 22.4));

    // IKs
    double axs = 4e-5 * (V - 18.80816) / (1.0 - exp(-(V - 18.80816) / 17.0));
    double bxs = 3.5e-5 * (V - 18.80816) / (exp((V - 18.80816) / 9.0) - 1.0);
    ALGEBRAIC[AV_iks_xs_inf] = 1.0 / sqrt(1.0 + exp(-(V - 19.9) / 12.7));
    ALGEBRAIC[AV_iks_xs_tau] = 0.5 / (axs + bxs);
    RATES[iks_xs] = (ALGEBRAIC[AV_iks_xs_inf] - STATES[iks_xs]) / ALGEBRAIC[AV_iks_xs_tau];
    ALGEBRAIC[AV_IKs] = CONSTANTS[AC_gKs] * pow(STATES[iks_xs], 2.0) * (V - EK);

    // IKur
    ALGEBRAIC[AV_gKur] = CONSTANTS[AC_gKur_amp] * (0.005 + 0.05 / (1.0 + exp(-(V - 15.0) / 13.0)));
    ALGEBRAIC[AV_ikur_ua_inf] = 1.0 / (1.0 + exp(-(V + 30.3) / 9.6));
    double aua = 0.65 / (exp(-(V + 10.0) / 8.5) + exp(-(V - 30.0) / 59.0));
    double bua = 0.65 / (2.5 + exp((V + 82.0) / 17.0));
    ALGEBRAIC[AV_ikur_ua_tau] = 1.0 / (CONSTANTS[AC_KQ10] * (aua + bua));
    RATES[ikur_ua] = (ALGEBRAIC[AV_ikur_ua_inf] - STATES[ikur_ua]) / ALGEBRAIC[AV_ikur_ua_tau];

    ALGEBRAIC[AV_ikur_uif_inf] = 1.0 / (1.0 + exp((V + 17.5) / 5.85));
    ALGEBRAIC[AV_ikur_uif_tau] = 400.0 + 1068.0 * exp(-pow(V / 50.0, 2.0));
    RATES[ikur_uif] = (ALGEBRAIC[AV_ikur_uif_inf] - STATES[ikur_uif]) / ALGEBRAIC[AV_ikur_uif_tau];

    ALGEBRAIC[AV_ikur_uis_inf] = ALGEBRAIC[AV_ikur_uif_inf];
    ALGEBRAIC[AV_ikur_uis_tau] = 2000.0 + 60000.0 * exp(-pow((V + 39.3) / 30.0, 2.0));
    RATES[ikur_uis] = (ALGEBRAIC[AV_ikur_uis_inf] - STATES[ikur_uis]) / ALGEBRAIC[AV_ikur_uis_tau];
    ALGEBRAIC[AV_IKur] = ALGEBRAIC[AV_gKur] * pow(STATES[ikur_ua], 3.0) * (0.25 * STATES[ikur_uif] + 0.75 * STATES[ikur_uis]) * (V - EK);

    // ICaL
    ALGEBRAIC[AV_ical_d_inf] = 1.0 / (1.0 + exp(-(V + 5.0) / 8.0));
    ALGEBRAIC[AV_ical_d_tau] = (fabs(V + 5.0) < 1e-6) ? 1.0 : (1.0 - exp(-(V + 5.0) / 6.24)) / (0.035 * (V + 5.0) * (1.0 + exp(-(V + 5.0) / 6.24)));
    RATES[ical_d] = (ALGEBRAIC[AV_ical_d_inf] - STATES[ical_d]) / ALGEBRAIC[AV_ical_d_tau];

    ALGEBRAIC[AV_ical_f_inf] = 1.0 / (1.0 + exp((V + 28.0) / 6.9));
    ALGEBRAIC[AV_ical_f_tau] = 9.0 / (0.0197 * exp(-0.03372 * pow(V + 5.0, 2.0)) + 0.02);
    RATES[ical_f] = (ALGEBRAIC[AV_ical_f_inf] - STATES[ical_f]) / ALGEBRAIC[AV_ical_f_tau];

    double fCa_inf = 1.0 / (1.0 + Cai / 0.00035);
    RATES[ical_fCa] = (fCa_inf - STATES[ical_fCa]) / CONSTANTS[AC_ical_fCa_tau];
    ALGEBRAIC[AV_ICaL] = CONSTANTS[AC_gCaL] * STATES[ical_d] * STATES[ical_f] * STATES[ical_fCa] * (V - 65.0);

    // INaK
    ALGEBRAIC[AV_fNaK] = 1.0 / (1.0 + 0.1245 * exp(-0.1 * V * CONSTANTS[AC_FRT]) + 0.0365 * CONSTANTS[AC_sigma] * exp(-V * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_INaK] = CONSTANTS[AC_INaK_max] * ALGEBRAIC[AV_fNaK] * CONSTANTS[AC_Ko] / (CONSTANTS[AC_Ko] + CONSTANTS[AC_KmKo]) / (1.0 + pow(CONSTANTS[AC_KmNai] / Nai, 1.5));

    // INaCa
    double B1_NaCa = exp(CONSTANTS[AC_gamma] * V * CONSTANTS[AC_FRT]) * pow(Nai, 3.0) * CONSTANTS[AC_Cao] - exp((CONSTANTS[AC_gamma] - 1.0) * V * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Nao], 3.0) * Cai;
    double B2_NaCa = (pow(CONSTANTS[AC_KmNa], 3.0) + pow(CONSTANTS[AC_Nao], 3.0)) * (CONSTANTS[AC_KmCa] + CONSTANTS[AC_Cao]) * (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_gamma] - 1.0) * V * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_INaCa] = CONSTANTS[AC_INaCa_max] * B1_NaCa / B2_NaCa;

    // Background Currents
    ALGEBRAIC[AV_IbNa] = CONSTANTS[AC_gbNa] * (V - ENa);
    ALGEBRAIC[AV_IbCa] = CONSTANTS[AC_gbCa] * (V - ECa);
    ALGEBRAIC[AV_IbK] = CONSTANTS[AC_gbK] * (V - EK);

    // IpCa
    ALGEBRAIC[AV_IpCa] = CONSTANTS[AC_IpCa_max] * Cai / (0.0005 + Cai);

    // Calcium Handling
    double Irel = CONSTANTS[AC_krel] * pow(STATES[ryr_u], 2.0) * STATES[ryr_v] * STATES[ryr_w] * (CaRel - Cai);
    double Itr = (CaUp - CaRel) / CONSTANTS[AC_tau_tr];
    double Iup = CONSTANTS[AC_I_up_max] / (1.0 + CONSTANTS[AC_K_up] / Cai);
    double Iup_leak = CONSTANTS[AC_I_up_max] * CaUp / CONSTANTS[AC_Ca_up_max];

    double Fn = 1e-12 * CONSTANTS[AC_V_rel] * Irel - 5e-13 / CONSTANTS[AC_F] * (0.5 * ALGEBRAIC[AV_ICaL] - 0.2 * ALGEBRAIC[AV_INaCa]) * CONSTANTS[AC_Cm];

    ALGEBRAIC[AV_iclca_qCa_inf] = 1.0 - pow(1.0 / (1.0 + Fn / 1.1e-10), 3.0);
    ALGEBRAIC[AV_iclca_qCa_tau] = 2.0;
    RATES[iclca_qCa] = (ALGEBRAIC[AV_iclca_qCa_inf] - STATES[iclca_qCa]) / ALGEBRAIC[AV_iclca_qCa_tau];
    ALGEBRAIC[AV_IClCa] = CONSTANTS[AC_gClCa] * STATES[iclca_qCa] * (V - ECl);

    ALGEBRAIC[AV_ryr_w_inf] = 1.0 - 1.0 / (1.0 + exp(-(V - 40.0) / 17.0));
    ALGEBRAIC[AV_ryr_w_tau] = 6.0 * (fabs(V - 7.9) < 1e-6 ? 2.0 / 13.0 : (1.0 - exp(-(V - 7.9) / 5.0)) / ((1.0 + 0.3 * exp(-(V - 7.9) / 5.0)) * (V - 7.9)));
    RATES[ryr_w] = (ALGEBRAIC[AV_ryr_w_inf] - STATES[ryr_w]) / ALGEBRAIC[AV_ryr_w_tau];

    double u_inf = 1.0 / (1.0 + exp(-(Fn - CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]));
    RATES[ryr_u] = (u_inf - STATES[ryr_u]) / CONSTANTS[AC_cajsr_u_tau];

    double v_inf = 1.0 - 1.0 / (1.0 + exp(-(Fn - 0.2 * CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]));
    double v_tau = 1.91 + 2.09 / (1.0 + exp(-(Fn - CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]));
    RATES[ryr_v] = (v_inf - STATES[ryr_v]) / v_tau;

    // dCa/dt
    double B1 = (2.0 * ALGEBRAIC[AV_INaCa] - (ALGEBRAIC[AV_IpCa] + ALGEBRAIC[AV_ICaL] + ALGEBRAIC[AV_IbCa])) * CONSTANTS[AC_Cm] / (2.0 * CONSTANTS[AC_V_i] * CONSTANTS[AC_F]) + (CONSTANTS[AC_V_up] * (Iup_leak - Iup) + Irel * CONSTANTS[AC_V_rel]) / CONSTANTS[AC_V_i];
    double B2 = 1.0 + CONSTANTS[AC_TRPN_max] * CONSTANTS[AC_Km_TRPN] / pow(Cai + CONSTANTS[AC_Km_TRPN], 2.0) + CONSTANTS[AC_CMDN_max] * CONSTANTS[AC_Km_CMDN] / pow(Cai + CONSTANTS[AC_Km_CMDN], 2.0);
    RATES[calcium_Cai] = B1 / B2;

    RATES[calcium_CaRel] = (Itr - Irel) * pow(1.0 + CONSTANTS[AC_CSQN_max] * CONSTANTS[AC_Km_CSQN] / pow(CaRel + CONSTANTS[AC_Km_CSQN], 2.0), -1.0);
    RATES[calcium_CaUp] = Iup - (Iup_leak + Itr * CONSTANTS[AC_V_rel] / CONSTANTS[AC_V_up]);

    // Concentrations d/dt
    RATES[sodium_Nai] = (-3.0 * ALGEBRAIC[AV_INaK] - (3.0 * ALGEBRAIC[AV_INaCa] + ALGEBRAIC[AV_IbNa] + ALGEBRAIC[AV_INa])) * CONSTANTS[AC_Cm] / (CONSTANTS[AC_V_i] * CONSTANTS[AC_F]);
    RATES[potassium_Ki] = (2.0 * ALGEBRAIC[AV_INaK] - (ALGEBRAIC[AV_IK1] + ALGEBRAIC[AV_IKur] + ALGEBRAIC[AV_IKr] + ALGEBRAIC[AV_IKs])) * CONSTANTS[AC_Cm] / (CONSTANTS[AC_V_i] * CONSTANTS[AC_F]);
    RATES[chloride_Cli] = ALGEBRAIC[AV_IClCa] * CONSTANTS[AC_Cm] / (CONSTANTS[AC_V_i] * CONSTANTS[AC_F]);

    // Total Iion
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_IK1] + ALGEBRAIC[AV_IKur] + ALGEBRAIC[AV_IKr] + ALGEBRAIC[AV_IKs] + ALGEBRAIC[AV_ICaL] + ALGEBRAIC[AV_IpCa] + ALGEBRAIC[AV_INaK] + ALGEBRAIC[AV_INaCa] + ALGEBRAIC[AV_IbNa] + ALGEBRAIC[AV_IbCa] + ALGEBRAIC[AV_IClCa];

    ALGEBRAIC[Istim] = PerisYague_2022computeIstim(VOI, stimulus);
    RATES[membrane_V] = 0.0;
    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }
}
