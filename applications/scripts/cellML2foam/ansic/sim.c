/*
Niederer_et_al_2006
Generated on 2026-01-30 15:19:07

Compiling on GCC:
 $ gcc -Wall -lm -lsundials_nvecserial -lsundials_cvode sim.c

Gnuplot example:
set terminal pngcairo enhanced linewidth 2 size 1200, 800;
set output 'V.png'
set size 1.0, 1.0
set xlabel 'time [ms]';
set grid
plot 'V.txt' using 1:2 with lines ls 1 title 'Vm'

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
  #include <sunmatrix/sunmatrix_dense.h>
  #include <sunlinsol/sunlinsol_dense.h>
  #include <cvodes/cvodes_direct.h>
#else
  #include <cvode/cvode_dense.h>
#endif

#define N_STATE 5

/* Declare intermediary, temporary and system variables */
static realtype t;
static realtype pace;
static realtype AC_A_1;
static realtype AC_A_2;
static realtype AC_A_3;
static realtype AV_Q;
static realtype AV_Tension;
static realtype AC_a;
static realtype AC_alpha_1;
static realtype AC_alpha_2;
static realtype AC_alpha_3;
static realtype AV_ExtensionRatio;
static realtype AC_dExtensionRatiodt;
static realtype AV_lambda;
static realtype AV_lambda_prev;
static realtype AV_time;
static realtype AC_beta_0;
static realtype AV_overlap;
static realtype AV_Ca_b;
static realtype AV_Ca_i;
static realtype AV_T_0;
static realtype AV_T_Base;
static realtype AC_T_ref;
static realtype AV_Ca_50;
static realtype AC_Ca_50ref;
static realtype AV_Ca_TRPN_50;
static realtype AC_K_1;
static realtype AC_K_2;
static realtype AC_K_z;
static realtype AC_alpha_0;
static realtype AV_alpha_Tm;
static realtype AC_alpha_r1;
static realtype AC_alpha_r2;
static realtype AC_beta_1;
static realtype AV_beta_Tm;
static realtype AC_n_Hill;
static realtype AC_n_Rel;
static realtype AV_z_max;
static realtype AC_z_p;
static realtype AC_Ca_TRPN_Max;
static realtype AV_J_TRPN;
static realtype AC_gamma_trpn;
static realtype AC_k_Ref_off;
static realtype AV_k_off;
static realtype AC_k_on;

/* Set values of constants */
static void
updateConstants(void)
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
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
    /* environment */
    ALGEBRAIC[AV_time] = t;
    
    /* Myofilaments */
    ALGEBRAIC[AV_ExtensionRatio] = ((ALGEBRAIC[AV_time] > 300000.0) ? 1.0 : 1.0);
    ALGEBRAIC[AV_lambda] = (((ALGEBRAIC[AV_ExtensionRatio] > 0.8) && (ALGEBRAIC[AV_ExtensionRatio] <= 1.15)) ? ALGEBRAIC[AV_ExtensionRatio] : ((ALGEBRAIC[AV_ExtensionRatio] > 1.15) ? 1.15 : 0.8));
    ALGEBRAIC[AV_lambda_prev] = ALGEBRAIC[AV_ExtensionRatio];
    
    /* filament_overlap */
    ALGEBRAIC[AV_overlap] = 1.0 + CONSTANTS[AC_beta_0] * (ALGEBRAIC[AV_lambda] - 1.0);
    
    /* Cross_Bridges */
    ALGEBRAIC[AV_Q] = STATES[2] + STATES[3] + STATES[4];
    RATES[2] = CONSTANTS[AC_A_1] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_1] * STATES[2];
    RATES[3] = CONSTANTS[AC_A_2] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_2] * STATES[3];
    RATES[4] = CONSTANTS[AC_A_3] * CONSTANTS[AC_dExtensionRatiodt] - CONSTANTS[AC_alpha_3] * STATES[4];
    
    /* intracellular_ion_concentrations */
    ALGEBRAIC[AV_Ca_i] = ((ALGEBRAIC[AV_time] < 1.0) ? 1000.0 * 1.84330000000000010e-07 : (((ALGEBRAIC[AV_time] >= 10.0) && (ALGEBRAIC[AV_time] < 15.0)) ? 1000.0 * (1.055 * pow(ALGEBRAIC[AV_time] / 1000.0, 3.0) - 0.03507 * pow(ALGEBRAIC[AV_time] / 1000.0, 2.0) + 0.0003992 * ALGEBRAIC[AV_time] / 1000.0 - 1.356e-06) : (((ALGEBRAIC[AV_time] >= 15.0) && (ALGEBRAIC[AV_time] < 55.0)) ? 1000.0 * (0.014 * pow(ALGEBRAIC[AV_time] / 1000.0, 3.0) - 0.002555 * pow(ALGEBRAIC[AV_time] / 1000.0, 2.0) + 0.0001494 * ALGEBRAIC[AV_time] / 1000.0 - 1.428e-06) : (((ALGEBRAIC[AV_time] >= 55.0) && (ALGEBRAIC[AV_time] < 250.0)) ? 1000.0 * (1.739e-05 * pow(ALGEBRAIC[AV_time] / 1000.0, 3.0) - 3.209e-06 * pow(ALGEBRAIC[AV_time] / 1000.0, 2.0) - 5.689e-06 * ALGEBRAIC[AV_time] / 1000.0 + 1.719e-06) : (((ALGEBRAIC[AV_time] >= 250.0) && (ALGEBRAIC[AV_time] < 490.0)) ? 1000.0 * (0.0001321 * pow(ALGEBRAIC[AV_time] / 1000.0, 4.0) - 0.0002197 * pow(ALGEBRAIC[AV_time] / 1000.0, 3.0) + 0.0001374 * pow(ALGEBRAIC[AV_time] / 1000.0, 2.0) - 3.895e-05 * ALGEBRAIC[AV_time] / 1000.0 + 4.441e-06) : 1000.0 * 1.21480000000000000e-07)))));
    
    /* tropomyosin */
    ALGEBRAIC[AV_Ca_50] = CONSTANTS[AC_Ca_50ref] * (1.0 + CONSTANTS[AC_beta_1] * (ALGEBRAIC[AV_lambda] - 1.0));
    ALGEBRAIC[AV_beta_Tm] = CONSTANTS[AC_alpha_r1] + CONSTANTS[AC_alpha_r2] * pow(STATES[1],
                                                                       CONSTANTS[AC_n_Rel] - 1.0) / (pow(STATES[1], CONSTANTS[AC_n_Rel]) + pow(CONSTANTS[AC_K_z], CONSTANTS[AC_n_Rel]));
    
    /* *remaining* */
    ALGEBRAIC[AV_Ca_b] = CONSTANTS[AC_Ca_TRPN_Max] - STATES[0];
    ALGEBRAIC[AV_Ca_TRPN_50] = ALGEBRAIC[AV_Ca_50] * CONSTANTS[AC_Ca_TRPN_Max] / (ALGEBRAIC[AV_Ca_50] + CONSTANTS[AC_k_Ref_off] / CONSTANTS[AC_k_on] * (1.0 - (1.0 + CONSTANTS[AC_beta_0] * (ALGEBRAIC[AV_lambda] - 1.0)) * 0.5 / CONSTANTS[AC_gamma_trpn]));
    ALGEBRAIC[AV_alpha_Tm] = CONSTANTS[AC_alpha_0] * pow(ALGEBRAIC[AV_Ca_b] / ALGEBRAIC[AV_Ca_TRPN_50],
                                                         CONSTANTS[AC_n_Hill]);
    ALGEBRAIC[AV_z_max] = (CONSTANTS[AC_alpha_0] / pow(ALGEBRAIC[AV_Ca_TRPN_50] / CONSTANTS[AC_Ca_TRPN_Max], CONSTANTS[AC_n_Hill]) - CONSTANTS[AC_K_2]) / (CONSTANTS[AC_alpha_r1] + CONSTANTS[AC_K_1] + CONSTANTS[AC_alpha_0] / pow(ALGEBRAIC[AV_Ca_TRPN_50] / CONSTANTS[AC_Ca_TRPN_Max], CONSTANTS[AC_n_Hill]));
    ALGEBRAIC[AV_T_Base] = CONSTANTS[AC_T_ref] * STATES[1] / ALGEBRAIC[AV_z_max];
    RATES[1] = ALGEBRAIC[AV_alpha_Tm] * (1.0 - STATES[1]) - ALGEBRAIC[AV_beta_Tm] * STATES[1];
    ALGEBRAIC[AV_T_0] = ALGEBRAIC[AV_T_Base] * ALGEBRAIC[AV_overlap];
    ALGEBRAIC[AV_Tension] = ((ALGEBRAIC[AV_Q] < 0.0) ? ALGEBRAIC[AV_T_0] * (CONSTANTS[AC_a] * ALGEBRAIC[AV_Q] + 1.0) / (1.0 - ALGEBRAIC[AV_Q]) : ALGEBRAIC[AV_T_0] * (1.0 + (CONSTANTS[AC_a] + 2.0) * ALGEBRAIC[AV_Q]) / (1.0 + ALGEBRAIC[AV_Q]));
    ALGEBRAIC[AV_k_off] = ((1.0 - ALGEBRAIC[AV_Tension] / (CONSTANTS[AC_gamma_trpn] * CONSTANTS[AC_T_ref]) > 0.1) ? CONSTANTS[AC_k_Ref_off] * (1.0 - ALGEBRAIC[AV_Tension] / (CONSTANTS[AC_gamma_trpn] * CONSTANTS[AC_T_ref])) : CONSTANTS[AC_k_Ref_off] * 0.1);
    ALGEBRAIC[AV_J_TRPN] = (CONSTANTS[AC_Ca_TRPN_Max] - STATES[0]) * ALGEBRAIC[AV_k_off] - ALGEBRAIC[AV_Ca_i] * STATES[0] * CONSTANTS[AC_k_on];
    RATES[0] = ALGEBRAIC[AV_J_TRPN];
    

    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    STATES[0] =  6.75931398649999987e-02;
    STATES[1] =  1.44179378369999993e-02;
    STATES[2] = 0.0;
    STATES[3] = 0.0;
    STATES[4] = 0.0;

}

/* Pacing event (non-zero stimulus) */
struct PacingEventS {
    double level;       /* The stimulus level (dimensionless, normal range [0,1]) */
    double start;       /* The time this stimulus starts */
    double duration;    /* The stimulus duration */
    double period;      /* The period with which it repeats (or 0 if it doesn't) */
    double multiplier;  /* The number of times this period occurs (or 0 if it doesn't) */
    struct PacingEventS* next;
};
typedef struct PacingEventS PacingEvent;

/*
 * Schedules a pacing event.
 * @param top The first event in a stack (the stack's head)
 * @param add The event to schedule
 * @return The new pointer to the head of the stack
 */
static PacingEvent*
PacingEvent_Schedule(PacingEvent* top, PacingEvent* add)
{
    add->next = 0;
    if (add == 0) return top;
    if (top == 0) return add;
    if (add->start <= top->start) {
        add->next = top;
        return add;
    }
    PacingEvent* evt = top;
    while(evt->next != 0 && evt->next->start <= add->start) {
        evt = evt->next;
    }
    add->next = evt->next;
    evt->next = add;
    return top;
}

/* CVODE Flags */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return(1);
    } /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
            return(1);
        }
    } /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return(1);
    }
    return 0;
}

/* Show output */
static void PrintOutput(realtype t, realtype y)
{
    #if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%4.1f     %14.6Le\n", t, y);
    #elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%4.1f     %14.6le\n", t, y);
    #else
        printf("%4.1f     %14.6e\n", t, y);
    #endif
    return;
}

/* Run a simulation */
int main()
{
    int success = -1;

    /* Sundials error flag */
    int flag;

    /* Declare variables that will need freeing */
    PacingEvent* events = NULL;
    N_Vector y = NULL;
    N_Vector dy = NULL;
    void *cvode_mem = NULL;

    #if SUNDIALS_VERSION_MAJOR >= 3
    SUNMatrix sundials_dense_matrix;
    SUNLinearSolver sundials_linear_solver;
    #endif

    #if SUNDIALS_VERSION_MAJOR >= 6
    /* Create sundials context */
    SUNContext sundials_context;
    flag = SUNContext_Create(NULL, &sundials_context);
    if (check_flag(&flag, "SUNContext_Create", 1)) goto error;

    /* Create state vectors */
    y = N_VNew_Serial(N_STATE, sundials_context);
    dy = N_VNew_Serial(N_STATE, sundials_context);
    #else
    /* Create state vectors */
    y = N_VNew_Serial(N_STATE);
    dy = N_VNew_Serial(N_STATE);
    #endif
    if (check_flag((void*)y, "N_VNew_Serial", 0)) goto error;
    if (check_flag((void*)dy, "N_VNew_Serial", 0)) goto error;

    /* Set calculated constants */
    updateConstants();

    /* Set initial values */
    default_initial_values(y);

    /* Set integration times */
    double tMin = 0;
    double tMax = 1000;
    double tLog = 0;

    /* Create pacing events */

    int nPacing = 1;
    events = (PacingEvent*)malloc(sizeof(PacingEvent)*nPacing);
    if (events == 0) goto error;
    int iPacing = 0;
    events[iPacing].level = 1.0;
    events[iPacing].start = 100.0;
    events[iPacing].duration = 0.5;
    events[iPacing].period = 1000.0;
    events[iPacing].multiplier = 0;
    iPacing++;

    /* Schedule events, make "next" point to the first event */
    PacingEvent* next = events;
    PacingEvent* fire = events + 1;
    for(iPacing=1; iPacing<nPacing; iPacing++) {
        next = PacingEvent_Schedule(next, fire++);
    }

    /* Fast forward events to starting time */
    double tNext = next->start;
    double tDown = 0.0;
    fire = 0;
    while (tNext <= tMin) {
        /* Event over? */
        if (fire != 0 && tNext >= tDown) {
            fire = 0;
        }
        /* New event? */
        if (next != 0 && tNext >= next->start) {
            fire = next;
            next = next->next;
            tDown = fire->start + fire->duration;
            if (fire->period > 0) {
                if (fire->multiplier != 1) {
                    if (fire->multiplier > 1) fire->multiplier--;
                    fire->start += fire->period;
                    next = PacingEvent_Schedule(next, fire);
                } else {
                    fire->period = 0;
                }
            }
        }
        /* Set next time */
        tNext = tMax;
        if (fire != 0 && tDown < tNext) tNext = tDown;
        if (next != 0 && next->start < tNext) tNext = next->start;
    }
    if (fire != 0) {
        pace = fire->level;
    } else {
        pace = 0.0;
    }

    /* Set simulation starting time */
    t = tMin;

    /* Create solver */
    #if SUNDIALS_VERSION_MAJOR >= 6
    cvode_mem = CVodeCreate(CV_BDF, sundials_context);
    #elif SUNDIALS_VERSION_MAJOR >= 4
    cvode_mem = CVodeCreate(CV_BDF);
    #else
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    #endif
    if (check_flag((void*)cvode_mem, "CVodeCreate", 0)) goto error;
    flag = CVodeInit(cvode_mem, rhs, t, y);
    if (check_flag(&flag, "CVodeInit", 1)) goto error;

    #if SUNDIALS_VERSION_MAJOR >= 6
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE, N_STATE, sundials_context);
        if (check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNLinSol_Dense(y, sundials_dense_matrix, sundials_context);
        if (check_flag((void *)sundials_linear_solver, "SUNLinSol_Dense", 0)) goto error;

        /* Attach the matrix and solver to cvode */
        flag = CVodeSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if (check_flag(&flag, "CVodeSetLinearSolver", 1)) goto error;
    #elif SUNDIALS_VERSION_MAJOR >= 4
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE, N_STATE);
        if (check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNLinSol_Dense(y, sundials_dense_matrix);
        if (check_flag((void *)sundials_linear_solver, "SUNLinSol_Dense", 0)) goto error;

        /* Attach the matrix and solver to cvode */
        flag = CVodeSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if (check_flag(&flag, "CVodeSetLinearSolver", 1)) goto error;
    #elif SUNDIALS_VERSION_MAJOR >= 3
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE,N_STATE);
        if(check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNDenseLinearSolver(y, sundials_dense_matrix);
        if(check_flag((void *)sundials_linear_solver, "SUNDenseLinearSolver", 0)) goto error;

        /* Attach the matrix and linear solver to cvode */
        flag = CVDlsSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) goto error;
    #else
        /* Create dense matrix for use in linear solves */
        flag = CVDense(cvode_mem, N_STATE);
        if (check_flag(&flag, "CVDense", 1)) goto error;
    #endif

    /* Set tolerances */
    double reltol = RCONST(1.0e-4);
    double abstol = RCONST(1.0e-6);
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) goto error;

    /* Go! */
    while(1) {
        if (tMax < tNext) tNext = tMax;
        flag = CVode(cvode_mem, tNext, y, &t, CV_ONE_STEP);
        if (check_flag(&flag, "CVode", 1)) break;
        if (flag == CV_SUCCESS) {
            /* Shot past next discontinuity? */
            if (t > tNext) {
                flag = CVodeGetDky(cvode_mem, tNext, 0, y);
                if (check_flag(&flag, "CVodeGetDky", 1)) goto error;
                t = tNext;
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
                /* Recalculate logging values at this point in time */
                rhs(t, y, dy, 0);
            }
            /* Event over? */
            if (fire != 0 && t >= tDown) {
                pace = 0;
                fire = 0;
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
            }
            /* New event? */
            if (next != 0 && t >= next->start) {
                fire = next;
                next = next->next;
                pace = fire->level;
                tDown = fire->start + fire->duration;
                if (fire->period > 0) {
                    if (fire->multiplier == 1) {
                        fire->period = 0;
                    } else {
                        if (fire->multiplier > 1) fire->multiplier--;
                        fire->start += fire->period;
                        next = PacingEvent_Schedule(next, fire);
                    }
                }
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
            }
            /* Set next time */
            tNext = tMax;
            if (fire != 0 && tDown < tNext) tNext = tDown;
            if (next != 0 && next->start < tNext) tNext = next->start;
            /* Log */
            if (t >= tLog) {
                /* Log current position */
                PrintOutput(t, STATES[0]);
            }
        }
        if (t >= tMax) break;
    }

    /* Success! */
    success = 0;

error:
    /* Free allocated space */
    free(events);

    /* Free CVODE space */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(dy);
    CVodeFree(&cvode_mem);
    #if SUNDIALS_VERSION_MAJOR >= 6
    SUNContext_Free(&sundials_context);
    #endif

    /* Return */
    return success;
}
