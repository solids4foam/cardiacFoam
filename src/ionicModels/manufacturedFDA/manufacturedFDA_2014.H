#include "manufacturedFDA_2014Names.H"

const char* manufacturedFDASTATES_NAMES[NUM_STATES] = {
    "V", //
    "u1", //
    "u2", //
    "u3"  //

};
const char* manufacturedFDAALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "Iion"
};


// Dim 1
// Dim 2
// Dim 3
void
manufacturedFDAinitConsts(double* CONSTANTS, double* RATES, double *STATES, int dimension)
{

    CONSTANTS[Cm]  = 2.0;
    CONSTANTS[Chi] = 3.0;

    if (dimension == 1)
    {
        CONSTANTS[Beta] = -1.1;
    }
    else if (dimension == 2)
    {
        CONSTANTS[Beta] = -5.9;
    }
    else if (dimension == 3)
    {
        CONSTANTS[Beta] = -8.6;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid mesh dimension detected!" << exit(Foam::FatalError);
    }

}



void
manufacturedFDAcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int dimension, bool solveVmWithinODESolver)
{
    // Local variables
    RATES[V] = 0.0;
    RATES[u1] = pow(STATES[u1] + STATES[u3] - STATES[V], 2.0) * pow(STATES[u2], 2.0) + 0.5 * (STATES[u1] + STATES[u3] - STATES[V]) * pow(STATES[u2], 2.0) * (STATES[V] - STATES[u3]);
    RATES[u2] = - (STATES[u1] + STATES[u3] - STATES[V]) *  pow(STATES[u2], 3);
    RATES[u3] = 0;
    ALGEBRAIC[Iion] =  - (CONSTANTS[Cm]/2.0) * (STATES[u1] + STATES[u3] - STATES[V])  * pow(STATES[u2], 2) * (STATES[V] - STATES[u3]) + (CONSTANTS[Beta] * (STATES[V] - STATES[u3])) / CONSTANTS[Chi];

    // std::cout
    //     << "u1 = " << u1
    //         << "\nu2 = " << u2
    //         << "\nu3 = " << u3
    //         << "\nV = " << V
    //         << "\nbeta = " << Beta
    //         << "\nChi2 = " << Chi
    //         << "\nCm = " << Cm << std::endl;
}
