/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    newManufacturedFoam

Description
    Solves the reaction-diffusion equation for muscle electrophysiology
    (monodomain), with a run-time selectable ionic model.

    This version is kept generic:
      - The solver only knows Vm is the transmembrane potential.
      - All ionic state indexing is encapsulated in the ionic model.
      - Manufactured-solution checking is done via a model-provided export
        function, not by indexing the state vector in the solver.

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tmanufacturedFDA.H"
#include "manufacturedSolutionHandler.H" // helper class for manufactured solution
#include "explicitLoopHandler.H"         // helper class for explicit Loop
#include "implicitLoopHandler.H"         // helper class for implicit Loop
#include "notsure.H"                     //test for implicit solver on ODE-PDE


#include "ionicModel.H"
#include "pimpleControl.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

//Initialization of Handler functions
manufacturedSolutionHandler msHandler(mesh, ionicModel());
explicitLoopHandler explicitHandler(mesh, ionicModel());
implicitLoopHandler implicitHandler(mesh, ionicModel());
notsureLoopHandler notsureHandler(mesh, ionicModel());

// External state storage: one N-state vector per cell.

const label nStates = ionicModel->nEqns();
Field<Field<scalar>> states
(
    mesh.nCells(),
    Field<scalar>(nStates, 0.0)
);
// -------------------------------------------------------------
// Mesh spacing (dx), dimension (dim), and MS-specific N
// -------------------------------------------------------------
scalar dx = Foam::cbrt(mesh.V().average().value());
int dim = 0;                                  // spatial dimension                           
// Determine mesh dimension once
{
    dim = mesh.nGeometricD();
    Info << "Mesh geometric dimension = " << dim << nl;
}

//Initialization
if (ionicModel->hasManufacturedSolution())
    msHandler.initializeManufactured(Vm, u1, u2, u3);
else
    msHandler.initializeNonManufactured();

dx  = msHandler.dx();
dim = msHandler.dim();

pimpleControl pimple(mesh);
const Switch solveExplicit(solutionVariablesMemory.lookup("solveExplicit"));
const scalar CFL = readScalar(solutionVariablesMemory.lookup("CFL"));

//test for different solver
Switch notsure(0);     // default = off
if (solutionVariablesMemory.found("notsure"))
{
    notsure = Switch(solutionVariablesMemory.lookup("notsure"));
}

// Global dt and nsteps
scalar dt = runTime.deltaTValue();
int nsteps = int( std::ceil(runTime.endTime().value() / dt));

Info << "Exporting fields:" << ionicModel->exportedFieldNames() << endl;



// =============== CASE 1: EXPLICIT SOLVER ====================

if (solveExplicit)
{
    explicitHandler.initializeExplicit
    (
        dt,
        nsteps,
        chi.value(),
        Cm.value(),
        conductivity,
        CFL,
        dx,
        dim
    );

    runTime.setDeltaT(dt);
    nsteps = int(std::ceil(runTime.endTime().value() / dt));

    while (runTime.loop())
    {
        const scalar t0 = runTime.value() - dt;

        explicitHandler.explicitLoop
        (
            t0,
            dt,
            Vm,
            Iion,
            states,
            externalStimulusCurrent,
            stimulusCellIDs,
            stimulusIntensity.value(),
            stimulusDuration.value(),
            chi,
            Cm,
            conductivity
        );
        ionicModel->exportStates(states, outFields);
        #include "updateActivationTimes.H"
        runTime.write();
    }
}

else if (notsure)   

// =============== CASE 2: NOT SURE SOLVER ====================
{
    Info << "\nUsing NOT-SURE solver mode\n" << endl;

    while (runTime.loop())
    {
        const scalar currentTime = runTime.value();
        const scalar deltaT      = runTime.deltaTValue();
        const scalar t0          = currentTime - deltaT;

        notsureHandler.notsureLoop
        (
            t0,
            deltaT,
            Vm,
            Iion,
            states,
            externalStimulusCurrent,
            stimulusCellIDs,
            stimulusIntensity.value(),
            stimulusDuration.value(),
            chi,
            Cm,
            conductivity,
            pimple
        );

        ionicModel->exportStates(states, outFields);
        #include "updateActivationTimes.H"
        runTime.write();
    }
}


// =============== CASE 3: IMPLICIT SOLVER ====================

else
{
    Info << "\nUsing implicit solver\n" << endl;

    while (runTime.loop())
    {
        const scalar currentTime = runTime.value();
        const scalar deltaT      = runTime.deltaTValue();
        const scalar t0          = currentTime - deltaT;

        implicitHandler.implicitLoop
        (
            t0,
            deltaT,
            Vm,
            Iion,
            states,
            externalStimulusCurrent,
            stimulusCellIDs,
            stimulusIntensity.value(),
            stimulusDuration.value(),
            chi,
            Cm,
            conductivity,
            pimple
        );

        ionicModel->exportStates(states, outFields);
        #include "updateActivationTimes.H"
        runTime.write();
    }
}

// Manufactured-solution post-processing
if (ionicModel->hasManufacturedSolution())
{
    msHandler.postProcess(Vm, u1, u2, u3, dt, nsteps, solveExplicit);
}

Info << "End" << nl << endl;
return 0;

}
