#include "heaviside.H"
#include "BuenoOrovio_2008Names.H"

const char* BuenoOrovioSTATES_NAMES[NUM_STATES] = {
    "u", // Membrane voltage
    "v", // Fast sodium activation variable
    "w", // Slow sodium activation variable
    "s"  // Slow inactivation variable

};
const char* BuenoOrovioALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "Jfi", // Fast sodium current
    "Jso", // Slow sodium current
    "Jsi",  // Slow inactivation current
    "Jion", // total current
    "Istim" // Stimulus current
};


// EPI tissueFlag = 1
// M tissueFlag = 2
// ENDO tissueFlag = 3
void
BuenoOrovioinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)
{

        CONSTANTS[uO] = 0;
        CONSTANTS[uU] = (tissueFlag == 1) ? 1.55 : (tissueFlag == 2) ? 1.61 : 1.56;
        CONSTANTS[thetaV] = 0.3;
        CONSTANTS[thetaW] = 0.13;
        CONSTANTS[thetaVMinus] = (tissueFlag == 1) ? 0.006 : (tissueFlag == 2) ? 0.1 : 0.2;
        CONSTANTS[thetaO] = (tissueFlag == 2) ? 0.005 : 0.006;

        CONSTANTS[tauV1Minus] = (tissueFlag == 1) ? 60 : (tissueFlag == 2) ? 80 : 75;
        CONSTANTS[tauV2Minus] = (tissueFlag == 1) ? 1150 : (tissueFlag == 2) ? 1.4506 : 10;
        CONSTANTS[tauVPlus] = 1.4506;

        CONSTANTS[tauW1Minus] = (tissueFlag == 1) ? 60 : (tissueFlag == 2) ? 70 : 6;
        CONSTANTS[tauW2Minus] = (tissueFlag == 1) ? 15 : (tissueFlag == 2) ? 8 : 140;

        CONSTANTS[kWMinus] = (tissueFlag == 1) ? 65 : 200;
        CONSTANTS[uWMinus] = (tissueFlag == 1) ? 0.03 : 0.016;

        CONSTANTS[tauWPlus] = (tissueFlag == 1) ? 200 : 280;

        CONSTANTS[tauFi] = (tissueFlag == 1) ? 0.11 : (tissueFlag == 2) ? 0.078 : 0.1;

        CONSTANTS[tauO1] = (tissueFlag == 1) ? 400 : (tissueFlag == 2) ? 410 : 470;
        CONSTANTS[tauO2] = (tissueFlag == 2) ? 7 : 6;

        CONSTANTS[tauSo1] = (tissueFlag == 1) ? 30.0181 : (tissueFlag == 2) ? 91 : 40;
        CONSTANTS[tauSo2] = (tissueFlag == 1) ? 0.9957 : (tissueFlag == 2) ? 0.8 : 1.2;
        CONSTANTS[kSo] = (tissueFlag == 1) ? 2.0458 : (tissueFlag == 2) ? 2.1 : 2;
        CONSTANTS[uSo] = (tissueFlag == 2) ? 0.6 : 0.65;

        CONSTANTS[tauS1] = 2.7342;
        CONSTANTS[tauS2] = (tissueFlag == 1) ? 16 : (tissueFlag == 2) ? 4 : 2;
        CONSTANTS[kS] = 2.0994;
        CONSTANTS[uS] = 0.9087;

        CONSTANTS[tauSi] = (tissueFlag == 1) ? 1.8875 : (tissueFlag == 2) ? 3.3849 : 2.9013;
        CONSTANTS[tauWInfty] = (tissueFlag == 1) ? 0.07 : (tissueFlag == 2) ? 0.01 : 0.0273;
        CONSTANTS[wInftyStar] = (tissueFlag == 1) ? 0.94 : (tissueFlag == 2) ? 0.5 : 0.78;

        // Initial values for the gating variables
        STATES[u] = 0;
        STATES[v] = 1.0;
        STATES[w] = 1.0;
        STATES[s] = 0.0;
};


void
BuenoOroviocomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
    using Foam::heaviside;
        //Time Constants currents
        double tauSo = CONSTANTS[tauSo1]
                    + (CONSTANTS[tauSo2] - CONSTANTS[tauSo1])
                    * (1.0 + tanh(CONSTANTS[kSo] * (STATES[u] - CONSTANTS[uSo]))) / 2.0;

        double tauO = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaO])) * CONSTANTS[tauO1]
                    + heaviside(STATES[u] - CONSTANTS[thetaO]) * CONSTANTS[tauO2];

        // Calculate the three currents
        ALGEBRAIC[Jfi] = -STATES[v] * heaviside(STATES[u] - CONSTANTS[thetaV])
                    * (STATES[u] - CONSTANTS[thetaV]) * (CONSTANTS[uU] - STATES[u]) / CONSTANTS[tauFi];

        ALGEBRAIC[Jso] = (STATES[u] - CONSTANTS[uO]) * (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW])) / tauO
                    + heaviside(STATES[u] - CONSTANTS[thetaW]) / tauSo;

        ALGEBRAIC[Jsi] = -heaviside(STATES[u] - CONSTANTS[thetaW])
                    * STATES[w] * STATES[s] / CONSTANTS[tauSi];

        ALGEBRAIC[Jion] = ALGEBRAIC[Jfi] + ALGEBRAIC[Jso] + ALGEBRAIC[Jsi];

        //std::cout << "CURRENTS " << __LINE__ << std::endl;
        // Time constants for gating variables
        double tauVMinus = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaVMinus])) * CONSTANTS[tauV1Minus]
                        + heaviside(STATES[u] - CONSTANTS[thetaVMinus]) * CONSTANTS[tauV2Minus];

        double tauWMinus = CONSTANTS[tauW1Minus]
                        + (CONSTANTS[tauW2Minus] - CONSTANTS[tauW1Minus])
                        * (1.0 + tanh(CONSTANTS[kWMinus] * (STATES[u] - CONSTANTS[uWMinus]))) / 2.0;

        double tauS = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW])) * CONSTANTS[tauS1]
                    + heaviside(STATES[u] - CONSTANTS[thetaW]) * CONSTANTS[tauS2];

        double vInfty = heaviside(CONSTANTS[thetaVMinus]- STATES[u]);

        double wInfty = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaO]))
                    * (1.0 - STATES[u] / CONSTANTS[tauWInfty])
                    + heaviside(STATES[u] - CONSTANTS[thetaO]) * CONSTANTS[wInftyStar];

        // Calculate the rates of change for the gating variables
        RATES[v] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaV]))
                * (vInfty - STATES[v]) / tauVMinus
                - heaviside(STATES[u] - CONSTANTS[thetaV]) * STATES[v] / CONSTANTS[tauVPlus];

        RATES[w] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW]))
                * (wInfty - STATES[w]) / tauWMinus
                - heaviside(STATES[u] - CONSTANTS[thetaW]) * STATES[w] / CONSTANTS[tauWPlus];

        RATES[s] = ((1.0 + tanh(CONSTANTS[kS] * (STATES[u] - CONSTANTS[uS]))) / 2.0 - STATES[s]) / tauS;

        RATES[u] = 0.0;



        //  S1–S2 STIMULUS LOGIC
        // -------------------------------------------------------------
        double I_stim = 0.0;

        double T1 = CONSTANTS[stim_period_S1];
        int    n1 = CONSTANTS[nstim1];

        double T2 = CONSTANTS[stim_period_S2];
        int    n2 = CONSTANTS[nstim2];

        double tStart = CONSTANTS[stim_start];
        double tDur   = CONSTANTS[stim_duration];
        double amp    = CONSTANTS[stim_amplitude] / 85.7;

        bool enable_S1 = (T1 > 0.0 && n1 > 0);

        if (enable_S1)
        {
            // Safe: only executed when T1 > 0 → no divide-by-zero
            double phase_S1 = VOI - floor(VOI / T1) * T1;

            bool in_S1 =
            (
                phase_S1 >= tStart &&
                phase_S1 <= tStart + tDur &&
                VOI <= tStart + T1 * n1
            );

            if (in_S1)
            {
                I_stim = amp;
            }
        }
        else
        {
            // S1 disabled: do nothing
        }

        bool enable_S2 = (T2 > 0.0 && n2 > 0);

        if (I_stim == 0.0 && enable_S2)
        {
            double tS1End = tStart + T1 * n1;   // safe even if T1 or n1 = 0

            bool in_S2_total_window =
            (
                VOI >= tS1End &&
                VOI <= tS1End + n2 * T2
            );

            if (in_S2_total_window)
            {
                // Safe: only executed when T2 > 0 → no divide-by-zero in fmod
                double phase_S2 = fmod(VOI - tS1End, T2);

                if (phase_S2 <= tDur)
                {
                    I_stim = amp;
                }
            }
        }

        // -------------------------------------------------------------
        ALGEBRAIC[Istim] = I_stim;


        if (solveVmWithinODESolver)
            {
                RATES[u] = -ALGEBRAIC[Jion] - ALGEBRAIC[Istim];
            }
}













