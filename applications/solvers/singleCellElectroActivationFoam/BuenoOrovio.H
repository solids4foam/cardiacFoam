#include "heaviside.H"

enum CONSTANTS_INDEX {
    uO,
    uU,
    thetaV,
    thetaW,
    thetaVMinus,
    thetaO,
    tauV1Minus,
    tauV2Minus,
    tauVPlus,
    tauW1Minus,
    tauW2Minus,
    kWMinus,
    uWMinus,
    tauWPlus,
    tauFi,
    tauO1,
    tauO2,
    tauSo1,
    tauSo2,
    kSo,
    uSo,
    tauS1,
    tauS2,
    kS,
    uS,
    tauSi,
    tauWInfty,
    wInftyStar,
    stimulusStart,        // Start time of pulse in ms
    stimulusDuration,      // Pulse duration in ms
    stimulusPeriod,    // Time between pulses in ms
    stimulusAmplitude, // Amplitude of stimulus
    NUM_CONSTANTS 
};
enum STATES_INDEX {
    u, // Membrane voltage
    v, // Fast sodium activation variable
    w, // Slow sodium activation variable
    s, // Slow inactivation variable
    NUM_STATES 
};
const char* STATES_NAMES[NUM_STATES] = {
    "u", // Membrane voltage
    "v", // Fast sodium activation variable
    "w", // Slow sodium activation variable
    "s"  // Slow inactivation variable
};
enum ALGEBRAIC_INDEX{
    Jfi, // Fast sodium current
    Jso, // Slow sodium current
    Jsi, // Slow inactivation current
    Istim, // Stimulus current
    NUM_ALGEBRAIC 
};
const char* ALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "Jfi", // Fast sodium current
    "Jso", // Slow sodium current
    "Jsi", // Slow inactivation current
    "Istim" // Stimulus current
};



// EPI tissueFlag = 1 
// M tissueFlag = 2
// ENDO tissueFlag = 3
void
BuenoOrovioinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag)
{
        //Stimulus
        CONSTANTS[stimulusStart]= 20.0;     // Start time of pulse in ms
        CONSTANTS[stimulusDuration]= 0.2;      // Pulse duration in ms
        CONSTANTS[stimulusPeriod]= 1000.0;    // Time between pulses in ms
        CONSTANTS[stimulusAmplitude]= -85;    // Amplitude of stimulus

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

        // Initial values for the hating variables
        STATES[u] = 0; 
        STATES[v] = 1.0;
        STATES[w] = 1.0;
        STATES[s] = 0.0;
};


void
BuenoOroviocomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag)
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


        // Time constants for gating variables
        double tauVMinus = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaVMinus])) * CONSTANTS[tauV1Minus]
                        + heaviside(STATES[u] - CONSTANTS[thetaVMinus]) * CONSTANTS[tauV2Minus];

        double tauWMinus = CONSTANTS[tauW1Minus]
                        + (CONSTANTS[tauW2Minus] - CONSTANTS[tauW1Minus])
                        * (1.0 + tanh(CONSTANTS[kWMinus] * (STATES[u] - CONSTANTS[uWMinus]))) / 2.0;

        double tauS = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW])) * CONSTANTS[tauS1]
                    + heaviside(STATES[u] - CONSTANTS[thetaW]) * CONSTANTS[tauS2];

        double vInfty = heaviside(CONSTANTS[thetaVMinus] - STATES[u]);
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


        ALGEBRAIC[Istim] = (
            (VOI - floor(VOI / CONSTANTS[stimulusPeriod]) * CONSTANTS[stimulusPeriod] >= CONSTANTS[stimulusStart]) &&
            (VOI - floor(VOI / CONSTANTS[stimulusPeriod]) * CONSTANTS[stimulusPeriod] <= CONSTANTS[stimulusStart] + CONSTANTS[stimulusDuration])
            ? CONSTANTS[stimulusAmplitude]
            : 0.0);
        // Dimensionless rate of Voltage u
        RATES[u] = -ALGEBRAIC[Jfi] - ALGEBRAIC[Jso] - ALGEBRAIC[Jsi] - ALGEBRAIC[Istim];

}












    