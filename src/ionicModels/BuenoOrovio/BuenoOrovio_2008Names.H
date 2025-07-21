#pragma once

enum CONSTANTS_INDEX{
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
enum STATES_INDEX{
    u, // Membrane voltage
    v, // Fast sodium activation variable
    w, // Slow sodium activation variable
    s, // Slow inactivation variable
    NUM_STATES 
};
enum ALGEBRAIC_INDEX{
    Jfi, // Fast sodium current
    Jso, // Slow sodium current
    Jsi, // Slow inactivation current
    Istim, // Stimulus current
    NUM_ALGEBRAIC 
};

void
BuenoOrovioinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag);
void
BuenoOroviocomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag);












    