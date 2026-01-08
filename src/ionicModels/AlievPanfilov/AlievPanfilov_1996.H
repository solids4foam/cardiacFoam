
#include "AlievPanfilov_1996Names.H"
#include "stimulusIO.H"

const char* AlievPanfilovSTATES_NAMES[NUM_STATES] = {
    "u", 
    "recovery_r"
};
const char* AlievPanfilovALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_eps",
    "Istim",
    "Iion_cm"
};


static inline const Foam::HashTable<Foam::wordList>&
AlievPanfilovDependencyMap()
{
    static const Foam::HashTable<Foam::wordList> dep
    (
        { 
        
        }
    );

    return dep;
}
inline Foam::scalar computeIstim(Foam::scalar t, const double* C)
{
    return Foam::stimulusIO::computeStimulus
    (
        t,
        C[stim_start],
        C[stim_period_S1],
        C[stim_duration],
        C[stim_amplitude],
        Foam::label(C[nstim1]),
        C[stim_period_S2],
        Foam::label(C[nstim2])
    );
}



void
AlievPanfilovinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)
{

    CONSTANTS[AC_Vp] = 20.0;
    CONSTANTS[AC_Vr] = -80.0;
    CONSTANTS[AC_Vth] = -70.0;
    CONSTANTS[AC_epsilon] = 0.01;
    CONSTANTS[AC_k] = 8.0;
    CONSTANTS[AC_mu1] = 0.2;
    CONSTANTS[AC_mu2] = 0.3;
    
    /* membrane_potential */
    CONSTANTS[AC_a] = (CONSTANTS[AC_Vth] - CONSTANTS[AC_Vr]) / (CONSTANTS[AC_Vp] - CONSTANTS[AC_Vr]);

    
    /* --- Initial values --- */
    STATES[u] = 0.05;
    STATES[recovery_r] = 0.0;
};


void
AlievPanfilovcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{

    /* recovery_variable */
    ALGEBRAIC[AV_eps] = CONSTANTS[AC_epsilon] + CONSTANTS[AC_mu1] * STATES[recovery_r] / (CONSTANTS[AC_mu2] + STATES[u]);
    RATES[recovery_r] = ALGEBRAIC[AV_eps] * (-STATES[recovery_r] - CONSTANTS[AC_k] * STATES[u] * (STATES[u] - (CONSTANTS[AC_a] + 1.0)));
    

    ALGEBRAIC[Iion_cm] = CONSTANTS[AC_k] * STATES[u] * (STATES[u] - CONSTANTS[AC_a]) * (STATES[u] - 1.0) + STATES[recovery_r] * STATES[u];

    //Tension Variables
    //STATES[active_tension_Ta] = 0;
    //CONSTANTS[AC_kTa] = 8.0;
    //CONSTANTS[AC_e0] = 1.0;
    //ALGEBRAIC[AV_e] = ((STATES[u] < 0.05) ? 1.0 * CONSTANTS[AC_e0] : 10.0 * CONSTANTS[AC_e0]);
    //RATES[active_tension_Ta] = ALGEBRAIC[AV_e] * (CONSTANTS[AC_kTa] * STATES[u] - STATES[active_tension_Ta]);

    
    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[u] = - ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}








