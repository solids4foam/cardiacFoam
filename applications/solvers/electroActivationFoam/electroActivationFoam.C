/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    electroActivationFoam

Description
 Generic model-agnostic solver for cardiac electrophysiology  based on the
 monodomain reaction–diffusion equation.

 The solver:

   • Treats Vm solely as the transmembrane potential; no assumptions are
     made about ionic state variables.

   • Delegates all ionic-model state indexing and ODE evaluation to the
     run-time selectable ionicModel.

   • Supports multiple time-integration strategies (explicit, implicit)
     via dedicated loop-handler classes.

   • Provides infrastructure for manufactured-solution verification through
     model-supplied export functions, avoiding any solver-side indexing of
     ionic states.

   • Stores ionic state vectors externally (one N-state vector per cell),
     ensuring complete separation between solver logic and model detail.

Authors
   Simão Nieto de Castro, UCD.
   Philip Cardiff, UCD.

 \*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tmanufacturedFDA.H"
#include "manufacturedSolutionHandler.H" // helper class for manufactured solution
#include "explicitLoopHandler.H"         // helper class for explicit Loop
#include "implicitLoopHandler.H"         // helper class for implicit Loop
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

    // Initialisation of handler functions
    manufacturedSolutionHandler msHandler(mesh, ionicModel());
    explicitLoopHandler explicitHandler(mesh, ionicModel());
    implicitLoopHandler implicitHandler(mesh, ionicModel());

    // External state storage: one N-state vector per cell.
    const label nStates = ionicModel->nEqns();
    Field<Field<scalar>> states
    (
        mesh.nCells(),Field<scalar>(nStates, 0.0)
    );

    // Initialisation
    pimpleControl pimple(mesh);

    // Structured mesh initialization
    scalar dx  = Foam::cbrt(mesh.V().average().value());
    int dim = mesh.nGeometricD();

    scalar dt = runTime.deltaTValue();
    int nsteps = int( std::ceil(runTime.endTime().value() / dt));

    if (ionicModel->hasManufacturedSolution())
    {
        msHandler.initializeManufactured(Vm, outFields, dx ,dim);
    }

    // Solution methodology flag
    const Switch solveExplicit(solutionVariablesMemory.lookup("solveExplicit"));

    // Extract the names of the fields to be exported
    const wordList exportNames = ionicModel->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }

    // =============== CASE 1: EXPLICIT SOLVER ====================
    if (solveExplicit)
    {
        const scalar CFL = readScalar(solutionVariablesMemory.lookup("CFL"));

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
        nsteps = int(std::ceil(runTime.endTime().value()/dt));

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
    // =============== CASE 2: IMPLICIT SOLVER ====================
    else
    {
        Info<< "\nUsing implicit solver\n" << endl;

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
        msHandler.postProcess(Vm, outFields, dt, nsteps, solveExplicit);
    }

    Info<< "End" << nl << endl;
    return 0;
}
