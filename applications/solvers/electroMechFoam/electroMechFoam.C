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

Solver
    electroMechFoam

Description
    Generic model-agnostic solver for cardiac electrophysiology
    based on the monodomain reaction–diffusion equation. The solver:

    - Delegates all ionic-model state indexing and ODE evaluation to the
      run-time selectable ionicModel.

    - Supports multiple time-integration strategies (explicit, implicit)
      via dedicated loop-handler classes.

    - Provides infrastructure for manufactured-solution verification through
      model-supplied export functions, avoiding any solver-side indexing of
      ionic states.

    - Stores ionic state vectors externally (one N-state vector per cell),
      ensuring complete separation between solver logic and model detail.

Authors
   Simão Nieto de Castro, UCD.
   Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"       
#include "implicitLoopHandler.H"         // helper class for implicit Loop
#include "explicitLoopHandler.H"         // helper class for explicit Loop
#include "ionicModel.H"
#include "activeTensionModel.H"
#include "activeTensionIO.H"
#include "pimpleControl.H"
#include "Field.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

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
    int nsteps = int(ceil(runTime.endTime().value()/dt));


    // Active tension arrays (default lambda = 1)
    scalarField lambda(mesh.nCells(), 1.0);
    scalarField Ta(mesh.nCells(), 0.0);

    // Extract the names of the fields to be exported
    const wordList exportNames = ionicModel->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }

    // Solution methodology flag
    const Switch solveExplicit(solutionVariablesMemory.lookup("solveExplicit"));

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
                stimulusCellIDsList,
                stimulusStartTimes,
                stimulusIntensity.value(),
                stimulusDuration.value(),
                chi,
                Cm,
                conductivity
            );

            ionicModel->exportStates(states, outFields);
            activeTensionModel->calculateTension
            (
                runTime.value(),
                dt,
                lambda,
                Ta
            );
            activeTensionIO::exportStateFields
            (
                Ta,
                activeTensionNames,
                activeTensionOutFields
            );

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
                stimulusCellIDsList,
                stimulusStartTimes,
                stimulusIntensity.value(),
                stimulusDuration.value(),
                chi,
                Cm,
                conductivity,
                pimple
            );

            ionicModel->exportStates(states, outFields);
            activeTensionModel->calculateTension
            (
                runTime.value(),
                deltaT,
                lambda,
                Ta
            );
            activeTensionIO::exportStateFields
            (
                Ta,
                activeTensionNames,
                activeTensionOutFields
            );

            #include "updateActivationTimes.H"

            runTime.write();
        }
    }

    runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}
