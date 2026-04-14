Here is the full Markdown document designed for your agent. You can save this directly as Cardiac_Sim_Templates.md.

Markdown
# Cardiac Electrophysiology OpenFOAM Agent Templates

This document contains the machine-readable lookup tables and OpenFOAM dictionary templates for cardiac electrophysiology simulations.
 
1. Agent knowlege on the ionicModel variables.  

# Dynamic Ionic Variable Extraction Protocol

*Agent Instructions: Do not hardcode export variables. When generating an `electroProperties` file for a specific `ionicModel`, you must dynamically query the simulation environment to find the valid variables.*

## Step 1: Run the Utility Tool
Before writing the `electroProperties` file, execute the built-in OpenFOAM utility via the shell:
```bash
it reads <model> from the dictionary with ionicModel <key>
listCellModelsVariables
and then the code extracts IONIC_SPECIFIC_EXPORTS

Target Template Block:

C++
    outputVariables 
    { 
        // Inject your parsed list into {{ IONIC_SPECIFIC_EXPORTS }}
        ionic { export (Vm {{ IONIC_SPECIFIC_EXPORTS }}); debug (Vm Iion); } 
    }
Result Example:

C++
    outputVariables 
    { 
        ionic { export (Vm Jsi Jfi Jso); debug (Vm Iion); } 
    }

### Why this makes your agent much smarter:
1. **Error Prevention:** It prevents the agent from guessing variables and causing OpenFOAM `FatalError` crashes.
2. **Future-Proof:** You can compile new ionic models into your C++ library, and the agent will instantly be able to simulate them without you changing its prompt.
3. **Self-Correction:** If the agent tries an invalid model name, `cellListModels` will likely throw an error listing all *available* models, which the agent can read and use to self-correct its mistake!

2. Universal System Dictionaries
Agent Instructions: Write these directly to the system/ directory. These contain the numerical solvers and schemes for all possible workflows. OpenFOAM will automatically select the necessary blocks.

system/fvSchemes

C++
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         backward; // 2nd order implicit for transient Vm
}
gradSchemes
{
    default         leastSquares;
}
divSchemes
{
    default         none;
    div(phiU,psi)   Gauss linear; // Used by Eikonal
}
laplacianSchemes
{
    default         Gauss linear corrected;
}
interpolationSchemes
{
    default         linear;
}
snGradSchemes
{
    default         corrected;
}
system/fvSolution

C++
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    "Vm|VmFinal|u|uFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-11;
        relTol          0.0;
    }
    "psi|psiFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-11;
        relTol          0.1;
    }
}

PIMPLE
{
    nOuterCorrectors    1000;
    residualControl
    {
        "Vm|psi"
        {
            relTol          0.0;
            tolerance       1e-6;
        }
    }
}

relaxationFactors
{
    equations
    {
        psi    0.9;
    }
}
3. Workflow-Specific electroProperties Templates
Agent Instructions: Inject one of the following templates into constant/electroProperties. Replace {{ PLACEHOLDERS }} dynamically.

Template A: 3D Monodomain + 1D Purkinje

C++
myocardiumSolver  monodomainSolver;

monodomainSolverCoeffs
{
    conductivity [-1 -3 3 0 0 2 0] (0.1334 0 0 0.01761 0 0.01761);
    chi          [0 -1 0 0 0 0 0]  140000;
    cm           [-1 -4 4 0 0 2 0] 0.01;
    
    ionicModel   {{ IONIC_MODEL_NAME }};
    tissue       epicardialCells;

    solutionAlgorithm  explicit;
    solver             RKF45;
    initialODEStep     1e-6;
    maxSteps           1000000000;
    
    electrophysicsAdvanceScheme  staggeredElectrophysicsAdvanceScheme;

    // Agent: Inject variables from JSON here, default to Vm
    outputVariables 
    { 
        ionic { export (Vm {{ IONIC_SPECIFIC_EXPORTS }}); debug (Vm Iion); } 
    }

    monodomainStimulus
    {
        stimulusLocationMin  (0 0 5.5e-3);
        stimulusLocationMax  (1.5e-3 1.5e-3 7e-3);
        stimulusDuration     [0 0 1 0 0 0 0]  2e-3;
        stimulusIntensity    [0 -3 0 0 0 1 0] 50000;
        stimulusStartTime    0.0;
    }

    conductionNetworkDomains
    {
        purkinjeNetwork
        {
            conductionSystemDomain  purkinjeNetworkModel;
            edges ( ( 0  1  0.003  3.0 ) ( 1  2  0.005  3.0 ) ( 1  3  0.005  3.0 ) );
            pvjNodes     ( 2  3 );
            pvjLocations ( (0.00245  0.00115  0.00025) (0.00745  0.00115  0.00025) );

            rootStimulus { startTime 0.0; duration 0.002; intensity 1000.0; }

            purkinjeNetworkModelCoeffs
            {
                conductionSystemSolver  monodomain1DSolver;
                ionicModel      {{ IONIC_MODEL_NAME }};
                tissue          endocardialCells;
                vm1DRest        -0.084;
                solver          RKF45;
                initialODEStep  1e-6;
                maxSteps        1000000000;
                chi             1400.0;
                cm              1.0e-2;
                outputVariables { export (Vm Icoupling); debug (Vm); }
            }
        }
    }

    domainCouplings
    {
        couplingA
        {
            electroDomainCoupler  reactionDiffusionPvjCoupler;
            conductionNetworkDomain purkinjeNetwork;
            rPvj                  500.0;
            pvjRadius             6e-4;
            couplingMode          bidirectional;
        }
    }
}
Template B: 3D Monodomain + Pseudo-ECG

C++
myocardiumSolver  monodomainSolver;

monodomainSolverCoeffs
{
    conductivity [-1 -3 3 0 0 2 0] (0.1334 0 0 0.01761 0 0.01761);
    chi          [0 -1 0 0 0 0 0]  140000;
    cm           [-1 -4 4 0 0 2 0] 0.01;
    
    ionicModel   {{ IONIC_MODEL_NAME }};
    tissue       epicardialCells;
    
    solutionAlgorithm  explicit;
    solver             RKF45;
    initialODEStep     1e-6;
    maxSteps           1000000000;
    electrophysicsAdvanceScheme  staggeredElectrophysicsAdvanceScheme;

    // Agent: Inject variables from JSON here, default to Vm
    outputVariables 
    { 
        ionic { export (Vm {{ IONIC_SPECIFIC_EXPORTS }}); debug (Vm Iion); } 
    }

    monodomainStimulus
    {
        stimulusLocationMin  (0 0 5.5e-3);
        stimulusLocationMax  (1.5e-3 1.5e-3 7e-3);
        stimulusDuration     [0 0 1 0 0 0 0]  2e-3;
        stimulusIntensity    [0 -3 0 0 0 1 0] 50000;
        stimulusStartTime    0.0;
    }
}

ecgDomains
{
    bodySurfaceECG
    {
        ecgSolver  pseudoECG;
        electrodePositions
        {
            lead_I   ( 0.1  0.0  0.0 );
            lead_II  ( 0.0  0.1  0.0 );
        }
    }
}
Template C: 3D Eikonal + 1D Eikonal

C++
myocardiumSolver  eikonalSolver;

eikonalSolverCoeffs
{
    conductivity  [-1 -3 3 0 0 2 0]  (0.1334 0 0 0.01761 0 0.01761);
    chi           [0 -1 0 0 0 0 0]   140000;
    cm            [-1 -4 4 0 0 2 0]  0.01;
    
    c0            [0 1 -1 0 0 0 0]   0.8;
    eikonalAdvectionDiffusionApproach  true;

    stimulusLocationMin  (0 0 5.5e-3);
    stimulusLocationMax  (1.5e-3 1.5e-3 7e-3);

    conductionNetworkDomains
    {
        purkinjeNetwork
        {
            conductionSystemDomain  purkinjeNetworkModel;
            edges ( ( 0  1  0.003  3.0 ) ( 1  2  0.005  3.0 ) ( 1  3  0.005  3.0 ) );
            pvjNodes     ( 2  3 );
            pvjLocations ( (0.00245  0.00115  0.00025) (0.00745  0.00115  0.00025) );
            
            rootStimulus { startTime 0.0; duration 0.002; intensity 1000.0; }

            purkinjeNetworkModelCoeffs
            {
                conductionSystemSolver  eikonalSolver;
                c0  [0 1 -1 0 0 0 0]    2.0;
                
                // Eikonal only uses psi (activation time), no ionic variables needed
                outputVariables { export (psi); debug (psi); }
            }
        }
    }

    domainCouplings
    {
        couplingA
        {
            electroDomainCoupler  eikonalPvjCoupler;
            conductionNetworkDomain purkinjeNetwork;
            pvjRadius             6e-4;
            couplingMode          unidirectional;
        }
    }
}
