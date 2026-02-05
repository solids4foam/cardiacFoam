#include "GoktepeKuhl_2004Names.H"

static const char* GoktepeKuhlSTATES_NAMES[NUM_STATES] = {
    "Ta"
};

static const char* GoktepeKuhlALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_e",
    "AV_Vm",
    "AV_u"
};


static inline const Foam::HashTable<Foam::wordList>&
GoktepeKuhlDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}
void
GoktepeKuhlinitConsts(double* CONSTANTS, double* RATES, double* STATES)
{

    CONSTANTS[AC_Vr] = -80.0;

    // Smooth damping tension parameters
    CONSTANTS[AC_eInfty] = 1.0;     // mV-1
    CONSTANTS[AC_e0] = 0.1;         // mV-1
    CONSTANTS[AC_eXi] = 1.0;        // mV-1
    CONSTANTS[AC_Vshift] = 0.0;

    CONSTANTS[AC_kTa] = 47.9;      // Kpa mV−1

   
    /* initial conditions */
    STATES[Ta] = 0.0;
    //STATES[AV_u]= 0.0;
}


void
GoktepeKuhlcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
    /* membrane_potential */
    
    //Tension Calculation
    ALGEBRAIC[AV_e] =
        CONSTANTS[AC_e0]
        + (CONSTANTS[AC_eInfty] - CONSTANTS[AC_e0])
        * exp(-exp(-CONSTANTS[AC_eXi] * (ALGEBRAIC[AV_Vm] - CONSTANTS[AC_Vshift])));

    RATES[Ta] = ALGEBRAIC[AV_e] * (CONSTANTS[AC_kTa] * ALGEBRAIC[AV_u] - STATES[Ta]);

    //ALGEBRAIC[AV_u]=0; // u is updated in the ionic model
    
    
 
}
