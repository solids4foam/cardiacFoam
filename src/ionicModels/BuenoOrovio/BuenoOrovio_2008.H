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

#include "heaviside.H"
#include "BuenoOrovio_2008Names.H"
#include "stimulusIO.H"

const char* BuenoOrovioSTATES_NAMES[NUM_STATES] = {
    "u", // Membrane voltage
    "v", // Fast sodium activation variable
    "w", // Slow sodium activation variable
    "s"  // Slow inactivation variable

};
const char* BuenoOrovioALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "Jfi",       // 0
    "Jso",       // 1
    "Jsi",       // 2
    "Jion",      // 3
    "Istim",     // 4

    "tauSo",     // 5
    "tauO",      // 6
    "tauVMinus", // 7
    "tauWMinus", // 8
    "tauS",      // 9
    "vInfty",    // 10
    "wInfty"     // 11
};

const char* BuenoOrovioCONSTANTS_NAMES[NUM_CONSTANTS] = {
    "uO",
    "uU",
    "thetaV",
    "thetaW",
    "thetaVMinus",
    "thetaO",
    "tauV1Minus",
    "tauV2Minus",
    "tauVPlus",
    "tauW1Minus",
    "tauW2Minus",
    "kWMinus",
    "uWMinus",
    "tauWPlus",
    "tauFi",
    "tauO1",
    "tauO2",
    "tauSo1",
    "tauSo2",
    "kSo",
    "uSo",
    "tauS1",
    "tauS2",
    "kS",
    "uS",
    "tauSi",
    "tauWInfty",
    "wInftyStar",
};

static inline const Foam::HashTable<Foam::wordList>&
BuenoOrovioDependencyMap()
{
    static const Foam::HashTable<Foam::wordList> dep
    (
        {
        // Canonical BO2008 names
        {"Jfi", Foam::wordList({"Jfi", "v", "RATES_v"})},
        {"Jso", Foam::wordList({"Jso", "tauO", "tauSo"})},
        {"Jsi", Foam::wordList({"Jsi", "w", "s", "RATES_s"})},
        {"Jion", Foam::wordList({"Jfi", "Jso", "Jsi", "Jion"})},
        {"v", Foam::wordList({"v", "vInfty", "tauVMinus", "RATES_v"})},
        {"w", Foam::wordList({"w", "wInfty", "tauWMinus", "RATES_w"})},
        {"s", Foam::wordList({"s", "tauS", "RATES_s"})},

        // CellML / alternative current naming fallbacks
        {"I_fi", Foam::wordList({"Jfi", "v", "RATES_v"})},
        {"I_so", Foam::wordList({"Jso", "tauO", "tauSo"})},
        {"I_si", Foam::wordList({"Jsi", "w", "s", "RATES_s"})},
        {"Iion", Foam::wordList({"Jfi", "Jso", "Jsi", "Jion"})},
        {"I_ion", Foam::wordList({"Jfi", "Jso", "Jsi", "Jion"})},
        {"Iion_cm", Foam::wordList({"Jfi", "Jso", "Jsi", "Jion"})}
        }
    );

    return dep;
}
inline Foam::scalar computeIstim
(
    Foam::scalar t,
    const Foam::StimulusProtocol& stimulus
)
{
    return Foam::stimulusIO::computeStimulus(t, stimulus);
}



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
BuenoOroviocomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
{
    using Foam::heaviside;
        //Time Constants currents
        ALGEBRAIC[tauSo] = CONSTANTS[tauSo1]
                    + (CONSTANTS[tauSo2] - CONSTANTS[tauSo1])
                    * (1.0 + tanh(CONSTANTS[kSo] * (STATES[u] - CONSTANTS[uSo]))) / 2.0;

        ALGEBRAIC[tauO] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaO])) * CONSTANTS[tauO1]
                    + heaviside(STATES[u] - CONSTANTS[thetaO]) * CONSTANTS[tauO2];

        // Calculate the three currents
        ALGEBRAIC[Jfi] = -STATES[v] * heaviside(STATES[u] - CONSTANTS[thetaV])
                    * (STATES[u] - CONSTANTS[thetaV]) * (CONSTANTS[uU] - STATES[u]) / CONSTANTS[tauFi];

        ALGEBRAIC[Jso] = (STATES[u] - CONSTANTS[uO]) * (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW])) / ALGEBRAIC[tauO]
                    + heaviside(STATES[u] - CONSTANTS[thetaW]) / ALGEBRAIC[tauSo];

        ALGEBRAIC[Jsi] = -heaviside(STATES[u] - CONSTANTS[thetaW])
                    * STATES[w] * STATES[s] / CONSTANTS[tauSi];

        ALGEBRAIC[Jion] = ALGEBRAIC[Jfi] + ALGEBRAIC[Jso] + ALGEBRAIC[Jsi];

        // Time constants for gating variables
        ALGEBRAIC[tauVMinus] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaVMinus])) * CONSTANTS[tauV1Minus]
                        + heaviside(STATES[u] - CONSTANTS[thetaVMinus]) * CONSTANTS[tauV2Minus];

        ALGEBRAIC[tauWMinus] = CONSTANTS[tauW1Minus]
                        + (CONSTANTS[tauW2Minus] - CONSTANTS[tauW1Minus])
                        * (1.0 + tanh(CONSTANTS[kWMinus] * (STATES[u] - CONSTANTS[uWMinus]))) / 2.0;

        ALGEBRAIC[tauS] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW])) * CONSTANTS[tauS1]
                    + heaviside(STATES[u] - CONSTANTS[thetaW]) * CONSTANTS[tauS2];

        ALGEBRAIC[vInfty] = heaviside(CONSTANTS[thetaVMinus]- STATES[u]);

        ALGEBRAIC[wInfty] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaO]))
                    * (1.0 - STATES[u] / CONSTANTS[tauWInfty])
                    + heaviside(STATES[u] - CONSTANTS[thetaO]) * CONSTANTS[wInftyStar];

        // Calculate the rates of change for the gating variables
        RATES[v] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaV]))
                * (ALGEBRAIC[vInfty] - STATES[v]) / ALGEBRAIC[tauVMinus]
                - heaviside(STATES[u] - CONSTANTS[thetaV]) * STATES[v] / CONSTANTS[tauVPlus];

        RATES[w] = (1.0 - heaviside(STATES[u] - CONSTANTS[thetaW]))
                * (ALGEBRAIC[wInfty] - STATES[w]) / ALGEBRAIC[tauWMinus]
                - heaviside(STATES[u] - CONSTANTS[thetaW]) * STATES[w] / CONSTANTS[tauWPlus];

        RATES[s] = ((1.0 + tanh(CONSTANTS[kS] * (STATES[u] - CONSTANTS[uS]))) / 2.0 - STATES[s]) / ALGEBRAIC[tauS];


        ALGEBRAIC[Istim] = computeIstim(VOI, stimulus);

        if (solveVmWithinODESolver)
            {
                RATES[u] = -ALGEBRAIC[Jion] - ALGEBRAIC[Istim];
            }
}



        //std::cout << "CURRENTS " << __LINE__ << std::endl;








