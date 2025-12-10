#pragma once

enum CONSTANTS_INDEX{
    Cm,
    Beta,
    Chi,
    NUM_CONSTANTS
};
enum STATES_INDEX{
    V,
    u1,
    u2,
    u3,
    NUM_STATES
};
enum ALGEBRAIC_INDEX{
    Iion,
    NUM_ALGEBRAIC
};

void
tmanufacturedFDAinitConsts(double* CONSTANTS, double* RATES, double *STATES, int dimension);
void
tmanufacturedFDAcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int dimension, bool solveVmWithinODESolver);












