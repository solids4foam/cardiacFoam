#include "dictionary.H"
#include "scalar.H"
#include "TNNP_2004Names.H"
#include "stimulusIO.H"


 static const char* TNNP_STATES_NAMES[17] = {
    "V",        // 0: membrane voltage (millivolt)
    "K_i",      // 1: intracellular potassium (millimolar)
    "Na_i",     // 2: intracellular sodium (millimolar)
    "Ca_i",     // 3: intracellular calcium (millimolar)
    "Xr1",      // 4: Xr1 gate (dimensionless)
    "Xr2",      // 5: Xr2 gate (dimensionless)
    "Xs",       // 6: Xs gate (dimensionless)
    "m",        // 7: m gate (dimensionless)
    "h",        // 8: h gate (dimensionless)
    "j",        // 9: j gate (dimensionless)
    "d",        // 10: d gate (dimensionless)
    "f",        // 11: f gate (dimensionless)
    "fCa",      // 12: fCa gate (dimensionless)
    "s",        // 13: s gate (dimensionless)
    "r",        // 14: r gate (dimensionless)
    "Ca_SR",    // 15: sarcoplasmic reticulum calcium (millimolar)
    "g"         // 16: g gate (dimensionless)
};
static const char* TNNP_ALGEBRAIC_NAMES[70] = {
    "Istim",                    // 0
    "xr1_inf",                  // 1
    "xr2_inf",                  // 2
    "xs_inf",                   // 3
    "m_inf",                    // 4
    "h_inf",                    // 5
    "j_inf",                    // 6
    "d_inf",                    // 7
    "f_inf",                    // 8
    "alpha_fCa",               // 9
    "s_inf",                    // 10
    "r_inf",                    // 11
    "g_inf",                    // 12
    "E_Na",                     // 13
    "alpha_xr1",               // 14
    "alpha_xr2",               // 15
    "alpha_xs",                // 16
    "alpha_m",                 // 17
    "alpha_h",                 // 18
    "alpha_j",                 // 19
    "alpha_d",                 // 20
    "tau_f",                    // 21
    "beta_fCa",                // 22
    "tau_s",                    // 23
    "tau_r",                    // 24
    "d_g",                      // 25
    "E_K",                      // 26
    "beta_xr1",                // 27
    "beta_xr2",                // 28
    "beta_xs",                 // 29
    "beta_m",                  // 30
    "beta_h",                  // 31
    "beta_j",                  // 32
    "gama_fCa",                // 33
    "beta_d",                  // 34
    "E_Ks",                     // 35
    "tau_xr1",                 // 36
    "tau_xr2",                 // 37
    "tau_xs",                  // 38
    "tau_m",                   // 39
    "tau_h",                   // 40
    "tau_j",                   // 41
    "gamma_d",                 // 42
    "fCa_inf",                 // 43
    "E_Ca",                     // 44
    "tau_d",                   // 45
    "d_fCa",                   // 46
    "alpha_K1",                // 47
    "beta_K1",                 // 48
    "xK1_inf",                 // 49
    "i_K1",                     // 50
    "i_Kr",                     // 51
    "i_Ks",                     // 52
    "i_Na",                     // 53
    "i_b_Na",                   // 54
    "i_CaL",                    // 55
    "i_b_Ca",                   // 56
    "i_to",                     // 57
    "i_NaK",                    // 58
    "i_NaCa",                   // 59
    "i_p_Ca",                   // 60
    "i_p_K",                    // 61
    "i_rel",                    // 62
    "i_up",                     // 63
    "i_leak",                   // 64
    "ddt_Ca_i_total",           // 65
    "ddt_Ca_sr_total",          // 66
    "f_JCa_i_free",             // 67
    "f_JCa_sr_free",            // 68
    "Iion_cm"                   // 69
};



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


static inline const Foam::HashTable<Foam::wordList>&
TNNPDependencyMap()
{
    static const Foam::HashTable<Foam::wordList> dep
    (
        {
            // Fast sodium current
            { "INa",
              Foam::wordList
              ({
                  "i_Na",
                  "m_inf","h_inf","j_inf",
                  "tau_m","tau_h","tau_j"
              }) },


            // L-type calcium current
            { "ICaL",
              Foam::wordList
              ({
                  "i_CaL",
                  "d_inf","f_inf","fCa_inf",
                  "tau_d","tau_f"
                  // tau_fCa is CONSTANT
              }) },

            // Transient outward potassium current
            { "Ito",
              Foam::wordList
              ({
                  "i_to",
                  "r_inf","s_inf",
                  "tau_r","tau_s"
              }) },

            
            // Rapid delayed rectifier potassium current
            { "IKr",
              Foam::wordList
              ({
                  "i_Kr",
                  "xr1_inf","xr2_inf",
                  "tau_xr1","tau_xr2"
                }) },

            // Slow delayed rectifier potassium current
            { "IKs",Foam::wordList ({"i_Ks","xs_inf","tau_xs"}) },

            // Inward rectifier potassium current
            { "IK1",Foam::wordList({"i_K1","xK1_inf"}) },

            // Sodium–potassium pump  
            { "INaK",Foam::wordList({"i_NaK"}) },

            
            // Sodium–calcium exchanger
            { "INaCa",Foam::wordList({"i_NaCa"}) },


            // Background currents
            { "IbNa",Foam::wordList({ "i_b_Na" }) },
            { "IbCa",Foam::wordList({ "i_b_Ca" }) },

            // Plateau / pump currents
            { "IpCa",Foam::wordList({ "i_p_Ca" }) },
            { "IpK", Foam::wordList({ "i_p_K" }) },

            // SR release gate (Ca-driven, not Vm sweep)
            { "RyR_g",Foam::wordList({"g_inf"}) },

            // ------------------------------------------------------------
            // Total ionic current
            // ------------------------------------------------------------
            { "Iion",
              Foam::wordList
              ({
                  "i_K1","i_to","i_Kr","i_Ks","i_CaL","i_Na","i_b_Na",
                  "i_NaK","i_NaCa","i_b_Ca", "i_p_K","i_p_Ca","Iion_cm"
              }) }
        }
    );

    return dep;
}



void
TNNPinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)

{
STATES[V] = -86.2;
CONSTANTS[R] = 8.314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96.485;
CONSTANTS[Cm] = 185; 
CONSTANTS[V_c] = 16404;


CONSTANTS[P_kna] = 0.03;
CONSTANTS[K_o] = 5.4;
CONSTANTS[Na_o] = 140;
STATES[K_i] = 138.3;
STATES[Na_i] = 11.6;
CONSTANTS[Ca_o] = 2;
STATES[Ca_i] = 0.0002;
CONSTANTS[g_K1] = 5.405;
CONSTANTS[g_Kr] = 0.096;
STATES[Xr1] = 0;
STATES[Xr2] = 1;
//Conditional tissue conductance Gks
CONSTANTS[g_Ks] = (tissueFlag == 2)
          ? 0.062
          : 0.245;

STATES[Xs] = 0;
CONSTANTS[g_Na] = 14.838;
STATES[m] = 0;
STATES[h] = 0.75;
STATES[j] = 0.75;
CONSTANTS[g_bNa] = 0.00029;
CONSTANTS[g_CaL] = 0.175;
STATES[d] = 0;
STATES[f] = 1;
STATES[fCa] = 1;
CONSTANTS[g_bCa] = 0.000592;
//Conditional tissue conductance Gto
CONSTANTS[g_to] = (tissueFlag == 3)
          ? 0.073
          : 0.294;

STATES[s] = 1;
STATES[r] = 0;
CONSTANTS[P_NaK] = 1.362;
CONSTANTS[K_mk] = 1;
CONSTANTS[K_mNa] = 40;
CONSTANTS[K_NaCa] = 1000;
CONSTANTS[K_sat] = 0.1;
CONSTANTS[alpha_NCX] = 2.5;
CONSTANTS[gamma_NCX] = 0.35;
CONSTANTS[Km_Ca] = 1.38;
CONSTANTS[Km_Nai] = 87.5;
CONSTANTS[g_pCa] = 0.825;
CONSTANTS[K_pCa] = 0.0005;
CONSTANTS[g_pK] = 0.0146;
STATES[Ca_SR] = 0.2;
STATES[g] = 1;
CONSTANTS[tau_g] = 2;
CONSTANTS[a_rel] = 0.016464;
CONSTANTS[b_rel] = 0.25;
CONSTANTS[c_rel] = 0.008232;
CONSTANTS[K_up] = 0.00025;
CONSTANTS[V_leak] = 8e-5;
CONSTANTS[Vmax_up] = 0.000425;
CONSTANTS[Buf_c] = 0.15;
CONSTANTS[K_buf_c] = 0.001;
CONSTANTS[Buf_sr] = 10;
CONSTANTS[K_buf_sr] = 0.3;
CONSTANTS[V_sr] = 1094;
CONSTANTS[tau_fCa] = 2.00000;

}

void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver)
{

//INa 
ALGEBRAIC[m_inf] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
RATES[m] = (ALGEBRAIC[m_inf] - STATES[m])/ALGEBRAIC[tau_m];
ALGEBRAIC[h_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[V])+ 310000.*exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
RATES[h] = (ALGEBRAIC[h_inf] - STATES[h])/ALGEBRAIC[tau_h];
ALGEBRAIC[j_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[V]) -  6.94800e-06*exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[V]))/(1.00000+exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[V]))/(1.00000+exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
RATES[j] = (ALGEBRAIC[j_inf] - STATES[j])/ALGEBRAIC[tau_j];

//ICaL
ALGEBRAIC[d_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
RATES[d] = (ALGEBRAIC[d_inf] - STATES[d])/ALGEBRAIC[tau_d];
ALGEBRAIC[f_inf] = 1.00000/(1.00000+exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*exp(- pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[V])/10.0000));
RATES[f] = (ALGEBRAIC[f_inf] - STATES[f])/ALGEBRAIC[tau_f];

ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
RATES[fCa] = (ALGEBRAIC[fCa_inf]>STATES[fCa]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_fCa]);

//Ito
//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1) 
    ? 1.00000 / (1.00000 + exp((STATES[V] + 20.0000) / 5.00000)) 
    : 1.10000 / (1.00000 + exp((STATES[V] + 28.0000) / 6.00000));
ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;
RATES[s] = (ALGEBRAIC[s_inf] - STATES[s])/ALGEBRAIC[tau_s];
ALGEBRAIC[r_inf] = 1.00000/(1.00000+exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*exp(- pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
RATES[r] = (ALGEBRAIC[r_inf] - STATES[r])/ALGEBRAIC[tau_r];

//IKr
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
RATES[Xr1] = (ALGEBRAIC[xr1_inf] - STATES[Xr1])/ALGEBRAIC[tau_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
RATES[Xr2] = (ALGEBRAIC[xr2_inf] - STATES[Xr2])/ALGEBRAIC[tau_xr2];

//IKs
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
RATES[Xs] = (ALGEBRAIC[xs_inf] - STATES[Xs])/ALGEBRAIC[tau_xs];

//SR calcium RyR gate g
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
RATES[g] = (ALGEBRAIC[g_inf]>STATES[g]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_g]);


ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (pow(CONSTANTS[Km_Nai], 3.00000)+pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
RATES[Na_i] = ( - (ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ 3.00000*ALGEBRAIC[i_NaK]+ 3.00000*ALGEBRAIC[i_NaCa])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);

ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

RATES[V] = 0.0;
if (solveVmWithinODESolver)
    {
        RATES[V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

RATES[K_i] = ( - ((ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_p_K]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[i_NaK])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);
ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*pow(STATES[Ca_SR], 2.00000))/(pow(CONSTANTS[b_rel], 2.00000)+pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+pow(CONSTANTS[K_up], 2.00000)/pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
RATES[Ca_i] =  ALGEBRAIC[ddt_Ca_i_total]*ALGEBRAIC[f_JCa_i_free];
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));
RATES[Ca_SR] =  ALGEBRAIC[ddt_Ca_sr_total]*ALGEBRAIC[f_JCa_sr_free];
}





void
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{

ALGEBRAIC[f_inf] = 1.00000/(1.00000+exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*exp(- pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+exp((25.0000 - STATES[V])/10.0000));

//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1) 
    ? 1.00000 / (1.00000 + exp((STATES[V] + 20.0000) / 5.00000)) 
    : 1.10000 / (1.00000 + exp((STATES[V] + 28.0000) / 6.00000));

ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*exp(- pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*exp(- pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;

ALGEBRAIC[r_inf] = 1.00000/(1.00000+exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*exp(- pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ pow((1.00000+exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
ALGEBRAIC[m_inf] = 1.00000/pow(1.00000+exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
ALGEBRAIC[h_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*exp( 0.0790000*STATES[V])+ 310000.*exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
ALGEBRAIC[j_inf] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*STATES[V]) -  6.94800e-06*exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*STATES[V]))/(1.00000+exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*exp( 0.0570000*STATES[V]))/(1.00000+exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
ALGEBRAIC[d_inf] = 1.00000/(1.00000+exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (pow(CONSTANTS[Km_Nai], 3.00000)+pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*pow(STATES[Ca_SR], 2.00000))/(pow(CONSTANTS[b_rel], 2.00000)+pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+pow(CONSTANTS[K_up], 2.00000)/pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);
}
