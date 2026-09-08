/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
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

Application
    pimpleHFDIBFoam.C

Group
    grpIncompressibleSolvers

Description
    A modified version of pimpleFoam using a modified version of the
    openHFDIBDEM immsersed boundary library.

Ported by Philip Cardiff
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "openHFDIBDEM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    IOdictionary HFDIBDEMDict
    (
        IOobject
        (
            "HFDIBDEMDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    #include "readDynMeshDict.H"

    openHFDIBDEM  HFDIBDEM(mesh);
    HFDIBDEM.initialize(lambda, U, refineF, maxRefinementLevel, runTime.timeName());

    #include "initialMeshRefinement.H"

    const scalar ibCouplingCoeff
    (
        HFDIBDEMDict.lookupOrDefault<scalar>("ibCouplingCoeff", 0.8)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        HFDIBDEM.createBodies(lambda, refineF);

        lambda.correctBoundaryConditions();
        f.correctBoundaryConditions();
        refineF.correctBoundaryConditions();
        surface.correctBoundaryConditions();

        HFDIBDEM.preUpdateBodies(lambda,f);

        lambda.correctBoundaryConditions();
        f.correctBoundaryConditions();
        refineF.correctBoundaryConditions();
        surface.correctBoundaryConditions();

        HFDIBDEM.recreateBodies(lambda, refineF);

        lambda.correctBoundaryConditions();
        f.correctBoundaryConditions();
        refineF.correctBoundaryConditions();
        surface.correctBoundaryConditions();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.solver.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    lambda *= 0.0;
                    HFDIBDEM.recreateBodies(lambda, refineF);
                }

                Ui.correctBoundaryConditions();
                lambda.correctBoundaryConditions();
                f.correctBoundaryConditions();
                refineF.correctBoundaryConditions();
                surface.correctBoundaryConditions();

                f *= lambda;
            }

            Ui.correctBoundaryConditions();
            lambda.correctBoundaryConditions();
            f.correctBoundaryConditions();
            refineF.correctBoundaryConditions();
            surface.correctBoundaryConditions();

            #include "UEqn.H"

            Ui.correctBoundaryConditions();
            lambda.correctBoundaryConditions();
            f.correctBoundaryConditions();
            refineF.correctBoundaryConditions();
            surface.correctBoundaryConditions();

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "updating HFDIBDEM" << endl;
        HFDIBDEM.postUpdateBodies(lambda, f);
        HFDIBDEM.addRemoveBodies(lambda,U,refineF);
        HFDIBDEM.updateDEM(lambda,refineF);
        Info << "updated HFDIBDEM" << endl;

        runTime.write();

        if (runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
