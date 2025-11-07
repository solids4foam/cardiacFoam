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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    electroActivationFoam

Description
    Solves the reaction-diffusion equation for muscle electrophysiology
    stemming from the mono-domain approach, where the ionic model is run-time
    selectable.

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ionicModelFDA.H"
//#include "fieldInit.H"
#include "pimpleControl.H"
//#include "manufacturedFields.H"
#include "manufacturedFDA.H"
//#include "testSolver.H"
#include "monodomain_fv_cell_centred.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));
    volScalarField Vm
    (
        IOobject
        (
            "Vm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Vm", dimVoltage, -0.084),
        "zeroGradient"
    );

     volScalarField u1
     (
         IOobject
         (
             "u1",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("u1", dimless, 0.0),
           "zeroGradient"
     );

     volScalarField u2
     (
         IOobject
         (
             "u2",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("u2", dimless, 0.0),
           "zeroGradient"
     );

     volScalarField u3
     (
         IOobject
         (
             "u3",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("u3", dimless, 0.0),
           "zeroGradient"
     );

// Create commonly used dimensionSets for convenience
// Read electoActivationProperties
Info<< "Reading electroActivationProperties\n" << endl;
IOdictionary electroActivationProperties
(
    IOobject
    (
        "electroActivationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

     volTensorField conductivity
     (
         IOobject
         (
             "conductivity",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE
         ),
         mesh,
         dimensionedTensor
         (
             "zero",
             pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
             tensor::zero
         )
     );

     // Check if we actually read it, otherwise assign constant value
     if (!conductivity.headerOk())
     {
         Info<< "conductivity not found on disk, using conductivity from "
             << electroActivationProperties.name() << nl << endl;

         conductivity =
             dimensionedTensor
             (
                 dimensionedSymmTensor
                 (
                     "conductivity",
                     pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                     electroActivationProperties
                 ) & tensor(I)
             );
     }
     else
     {
         Info<< "conductivity field read from " << runTime.timeName() << nl << endl;
     }

    volScalarField Iion
    (
        IOobject
        (
            "ionicCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
        "zeroGradient"
    );

    // Read ionic model dict
    const dictionary& ionicModelCoeffs =
        electroActivationProperties.subDict("ionicModelCoeffs");

    // Create ionicModelCellML object
    autoPtr<ionicModelFDA> ionicModelFDA =
        ionicModelFDA::New(ionicModelCoeffs, mesh.nCells(), runTime.deltaTValue());


    // Initialise fields
    ionicModelFDA->initializeFields
    (
        Vm,
        u1,
        u2,
        u3,
        mesh.C()
    );

    const Switch solveExplicit
    (
        electroActivationProperties.lookup("solveExplicit")
    );

    //std::vector<int> Ns = {50, 100, 200, 400};

    //std::vector<RunResult> results;
    // results.reserve(Ns.size());
    // for (int N : Ns) results.push_back(run_solver_ccfv(N));

    //RunResult r = run_solver_ccfv(mesh.nCells());

    //
    // run solver function
    //
    
    const int N = mesh.nCells();
    const scalar cfl = 0.1;

    // Cell-centred grid on [0,1]: centres x_i = (i+0.5)*dx
    // const double dx = 1.0 / double(N);
    // Foam::scalarField x(N);
    // for (int i=0; i<N; ++i) x[i] = (i + 0.5)*dx;
    const scalarField x(mesh.C().component(vector::X));
    const double dx = 1.0 / double(N);

    // Initial conditions
    //Foam::scalarField V(N), u1(N), u2(N), u3(N, 0.0);
    for (int i=0; i<N; ++i)
    {
        Vm[i] = manufactured_V(0.0, x[i]);
        double a,b,c;
        manufactured_u(0.0, x[i], a, b, c);
        u1[i] = a; u2[i] = b; u3[i] = c;
    }
    Vm.correctBoundaryConditions();
    u1.correctBoundaryConditions();
    u2.correctBoundaryConditions();
    u3.correctBoundaryConditions();

    Vm.storeOldTime();
    u1.storeOldTime();
    u2.storeOldTime();
    u3.storeOldTime();

    const dimensionedScalar CmDS("Cm", dimCurrent*dimTime/(dimVoltage*dimArea), Foam::Cm);
    const dimensionedScalar chiDS("chi", dimArea/dimVolume, chi);
    //const dimensionedScalar sigmaDS("sigma", pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume), sigma);

    scalar dt = runTime.deltaTValue();
    int nsteps = int(std::ceil(runTime.endTime().value() / dt));

    if (solveExplicit)
    {
        // Explicit diffusion stability
        dt = cfl * dx*dx / max(conductivity.component(tensor::XX)).value(); // sigma;
        // int nsteps = int(std::ceil(Tfinal / dt));
        // dt = Tfinal / double(nsteps); // snap to hit Tfinal
        nsteps = int(std::ceil(runTime.endTime().value() / dt));
        dt = runTime.endTime().value() / double(nsteps); // snap to hit Tfinal

        runTime.setDeltaT(dt);

        // Time loop: forward Euler
        // Foam::scalarField L(N);
        // for (int step=0; step<nsteps; ++step)
        while (runTime.loop())
        {
            // Info<< "Time = " << runTime.value() << endl;
        //     runTime++;

            // Laplacian(V)
            //laplacian_neumann(Vm, dx, L);
            // L = fvc::laplacian(conductivity, Vm);

            // chi*Cm*dV/dt = sigma*L - chi*Iion
            // dV/dt = (sigma*L - chi*Iion)/(chi*Cm)
            // Foam::scalarField dVdt(N), Iion(N);
            // Foam::scalarField dVdt(N);
            //Foam::scalarField Iion(N);
            for (int i=0; i<N; ++i)
            {
                Iion[i] = I_ion(u1[i], u2[i], u3[i], Vm[i]);
                // // dVdt[i] = (sigma*L[i] - chi*Iion[i]) / (chi*Cm);
                // dVdt[i] = (L[i] - chi*Iion[i]) / (chi*Cm);
            }
            Iion.correctBoundaryConditions();

            // Update V explicitly
            // Foam::scalarField Vcorrect(Vm);
            // for (int i=0; i<N; ++i) Vcorrect[i] += dt * dVdt[i];
            // Vm.correctBoundaryConditions();
            //Info<< "Vcorrect = " << Vcorrect << endl;

            // Use OF eqn
            solve
            (
                chiDS*CmDS*fvm::ddt(Vm)
             == fvc::laplacian(conductivity, Vm)
              - chiDS*Iion
            );

            // Update u explicitly using the (new) V (matches the Python ordering)
            for (int i=0; i<N; ++i)
            {
                double rhs1, rhs2, rhs3;
                f_rhs(u1[i], u2[i], u3[i], Vm[i], rhs1, rhs2, rhs3);
                u1[i] += dt * rhs1;
                u2[i] += dt * rhs2;
                // u3 stays 0
            }
            u1.correctBoundaryConditions();
            u2.correctBoundaryConditions();
            //u3.correctBoundaryConditions();
        }
    }
    else // solve implicit
    {
        // Info<< "deltaT (before) = " << runTime.deltaTValue() << endl;
        // runTime.setDeltaT(dt);
        // Info<< "deltaT (after) = " << runTime.deltaTValue() << endl;

        while (runTime.loop())
        {
            Info<< nl << "Time = " << runTime.value() << endl;

            // Update the old-time STATES within the ionic model
            refCast<manufacturedFDA>(*ionicModelFDA).updateStatesOld();

            // Update ionic current explicitly
            for (int i=0; i<N; ++i)
            {
                Iion[i] = I_ion(u1[i], u2[i], u3[i], Vm[i]);
            }
            Iion.correctBoundaryConditions();

            // Outer iteraiton implicit loop
            while (pimple.loop())
            {
                // Update Vm
                solve
                (
                    chiDS*CmDS*fvm::ddt(Vm)
                 == fvm::laplacian(conductivity, Vm)
                  - chiDS*Iion
                );
            }

            // Update ionic model explicitly

            // Update the old-time STATES within the ionic model
            //// Before solving the ionic model, reset its state to the old time
            //refCast<manufacturedFDA>(*ionicModelFDA).resetStatesToStatesOld();
            ionicModelFDA->calculateCurrent
            (
                runTime.value() - runTime.deltaTValue(),
                runTime.deltaTValue(),
                Vm.internalField(),
                Iion,
                u1,
                u2,
                u3
            );
            Vm.correctBoundaryConditions();
            u1.correctBoundaryConditions();
            u2.correctBoundaryConditions();
            u3.correctBoundaryConditions();

            // // Update u explicitly using the (new) V (matches the Python ordering)
            // // We could use sub-stepping here - we don't need to use the same
            // // dt as the PDE above => actually the ODE class already does this!
            // for (int i=0; i<N; ++i)
            // {
            //     double rhs1, rhs2, rhs3;
            //     f_rhs(u1[i], u2[i], u3[i], Vm[i], rhs1, rhs2, rhs3);
            //     // u1[i] += dt * rhs1;
            //     // u2[i] += dt * rhs2;
            //     u1[i] = u1.oldTime()[i] + dt*rhs1;
            //     u2[i] = u2.oldTime()[i] + dt*rhs2;
            //     // u3 stays 0
            // }
            // u1.correctBoundaryConditions();
            // u2.correctBoundaryConditions();
            //u3.correctBoundaryConditions();
        }
    }

    // Errors vs manufactured solution at Tfinal
    Foam::scalarField Vex(N), u1e(N), u2e(N), u3e(N);
    for (int i=0; i<N; ++i)
    {
        Vex[i] = manufactured_V(Tfinal, x[i]);
        double a,b,c;
        manufactured_u(Tfinal, x[i], a, b, c);
        u1e[i]=a; u2e[i]=b; u3e[i]=c;
    }
    auto eV  = error_norms(Vm,   Vex);
    auto eU1 = error_norms(u1,  u1e);
    auto eU2 = error_norms(u2,  u2e);

    RunResult r;
    r.N = N; r.dx = dx; r.dt = dt; r.steps = nsteps;
    r.errV_inf  = eV.Linf;
    r.errU1_inf = eU1.Linf;
    r.errU2_inf = eU2.Linf;






    std::cout.setf(std::ios::fixed);
    std::cout << "   N      dx        dt       steps      errV_inf      pV"
              << "      errU1_inf     pU1      errU2_inf     pU2\n";

    // for (size_t i=0; i<results.size(); ++i)
    {
        //const auto& r = results[i];
        double pV  = NAN, pU1 = NAN, pU2 = NAN;
        // const int i = 0;
        // if (i>0)
        // {
            // //const auto& rc = results[i-1];
            // const auto& rc = r;
            // pV  = observed_order(rc.errV_inf,  r.errV_inf,  rc.dx, r.dx);
            // pU1 = observed_order(rc.errU1_inf, r.errU1_inf, rc.dx, r.dx);
            // pU2 = observed_order(rc.errU2_inf, r.errU2_inf, rc.dx, r.dx);
        // }

        std::cout << std::setw(4) << r.N
                  << std::setw(9) << std::setprecision(6) << r.dx
                  << std::setw(11) << std::setprecision(6) << r.dt
                  << std::setw(10) << r.steps
                  << std::setw(14) << std::setprecision(9) << r.errV_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pV)?0.0:pV)
                  << std::setw(14) << std::setprecision(9) << r.errU1_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pU1)?0.0:pU1)
                  << std::setw(14) << std::setprecision(9) << r.errU2_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pU2)?0.0:pU2)
                  << "\n";
    }

    // std::cout << "\nNotes:\n"
    //           << " - Cell-centred FV grid: x_i = (i+0.5)*dx on [0,1].\n"
    //           << " - Explicit Euler in time with dt ~ dx^2 for diffusion stability.\n";


    // // pimpleControl pimple(mesh);

    // // #include "createFields.H"
    // // Transmembrane potential

    // computeAndPrintErrors
    // (
    //     Vm, u1, u2, u3, mesh.C().component(vector::X), runTime.value()
    // );


    // computeAndPrintErrors
    // (
    //     Vm, u1, u2, u3, mesh.C().component(vector::X), runTime.value()
    // );

    // // Initial write of fields
    // // Force write of initialized fields at t=0
    // Vm.write();
    // u1.write();
    // u2.write();
    // u3.write();
    // runTime.write();
    // Info<< "conductivity[0] = " << conductivity[0] << nl
    // Info<< "conductivity = " << conductivity << nl
    //     << "chi = " << chi << nl
    //     << "Cm = " << Cm << endl;

    // Info<< "Initial fields written to 0/ directory\n" << endl;
    // {
    //     double pV  = NAN, pU1 = NAN, pU2 = NAN;
    //     if (i>0)
    //     {
    //         const auto& rc = results[i-1];
    //         pV  = observed_order(rc.errV_inf,  r.errV_inf,  rc.dx, r.dx);
    //         pU1 = observed_order(rc.errU1_inf, r.errU1_inf, rc.dx, r.dx);
    //         pU2 = observed_order(rc.errU2_inf, r.errU2_inf, rc.dx, r.dx);
    //     }

    //     std::cout << std::setw(4) << r.N
    //               << std::setw(9) << std::setprecision(6) << r.dx
    //               << std::setw(11) << std::setprecision(6) << r.dt
    //               << std::setw(10) << r.steps
    //               << std::setw(14) << std::setprecision(6) << r.errV_inf
    //               << std::setw(9)  << std::setprecision(3) << (std::isnan(pV)?0.0:pV)
    //               << std::setw(14) << std::setprecision(6) << r.errU1_inf
    //               << std::setw(9)  << std::setprecision(3) << (std::isnan(pU1)?0.0:pU1)
    //               << std::setw(14) << std::setprecision(6) << r.errU2_inf
    //               << std::setw(9)  << std::setprecision(3) << (std::isnan(pU2)?0.0:pU2)
    //               << "\n";
    // }


    // // Cell-centred grid on [0,1]: centres x_i = (i+0.5)*dx
    // const label N = mesh.nCells();
    // // const double dx = 1.0 / double(N);
    // // std::vector<double> x(N);
    // // for (int i=0; i<N; ++i) x[i] = (i + 0.5)*dx;

    // // Explicit diffusion stability
    // // const double dt = cfl * dx*dx / sigma;
    // scalar dt = runTime.deltaTValue();
    // int nsteps = int(std::ceil(runTime.endTime().value() / dt));
    // //dt = Tfinal / double(nsteps); // snap to hit Tfinal

    // // Cell centre coordinates
    // const scalarField x(mesh.C().component(Foam::vector::X));

    // // Initial conditions
    // scalarField V(N), u1(N), u2(N), u3(N, 0.0);
    // //scalarField u1(N), u2(N), u3(N, 0.0);
    // for (int i=0; i<N; ++i)
    // {
    //     Vm[i] = manufactured_V(0.0, x[i]);
    //     scalar a,b,c;
    //     manufactured_u(0.0, x[i], a, b, c);
    //     u1[i] = a; u2[i] = b; u3[i] = c;
    // }
    // // Vm.primitiveFieldRef() = V;
    // Vm.correctBoundaryConditions();

    // const dimensionedScalar CmDS("Cm", dimCurrent*dimTime/(dimVoltage*dimArea), Cm);
    // const dimensionedScalar chiDS("chi", dimArea/dimVolume, chi);
    // const dimensionedScalar sigmaDS("sigma", pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume), sigma);

    // // Loop through time
    // Info<< "\nStarting time loop\n" << endl;

    // scalarField L(N);
    // for (int step=0; step<nsteps; ++step)
    // {
    //     // Laplacian(V)
    //     laplacian_neumann(V, dx, L);

    //     // dV/dt = (sigma*L - chi*Iion)/(chi*Cm)
    //     scalarField dVdt(N), Iion(N);
    //     for (int i=0; i<N; ++i)
    //     {
    //         Iion[i] = I_ion(u1[i], u2[i], u3[i], Vm[i]);
    //         dVdt[i] = (sigma*L[i] - chi*Iion[i]) / (chi*Cm);
    //     }

    //     // Update V explicitly
    //     // for (int i=0; i<N; ++i) V[i] += dt * dVdt[i];
    //     Vm.primitiveFieldRef() += dt*dVdt;
    //     Vm.correctBoundaryConditions();
    //     // Update u explicitly using the (new) V (matches the Python ordering)

    //     for (int i=0; i<N; ++i)
    //     {
    //         double rhs1, rhs2, rhs3;
    //         f_rhs(u1[i], u2[i], u3[i], Vm[i], rhs1, rhs2, rhs3);
    //         u1[i] += dt * rhs1;
    //         u2[i] += dt * rhs2;
    //         // u3 stays 0
    //     }
    // }

    // // while (runTime.loop())
    // // {
    //     //Info<< nl << "Time = " << runTime.timeName() << endl;

    //     // OF Vm update

    //     // forAll(ionicCurrent, cellI)
    //     // {
    //     //     ionicCurrent[cellI] =
    //     //         I_ion(u1[cellI], u2[cellI], u3[cellI], Vm[cellI]);
    //     // }
    //     // ionicCurrent.correctBoundaryConditions();

    //     // fvScalarMatrix VmEqn
    //     // (
    //     //     fvm::ddt(Vm)
    //     //  == (sigmaDS/(chiDS*CmDS))*fvm::laplacian(Vm) - ionicCurrent/CmDS
    //     //  // == (sigmaDS/(chiDS*CmDS))*fvc::laplacian(Vm) - ionicCurrent/CmDS
    //     // );
    //     // VmEqn.solve();

    // //     // Update u explicitly using the (new) V (matches the Python ordering)
    // //     for (int i=0; i<N; ++i)
    // //     {
    // //         double rhs1, rhs2, rhs3;
    // //         f_rhs(u1[i], u2[i], u3[i], Vm[i], rhs1, rhs2, rhs3);
    // //         u1[i] += dt * rhs1;
    // //         u2[i] += dt * rhs2;
    // //         // u3 stays 0
    // //     }
    // // }

    // // Errors vs manufactured solution at Tfinal
    // scalarField Vex(N), u1e(N), u2e(N), u3e(N);
    // for (int i=0; i<N; ++i)
    // {
    //     Vex[i] = manufactured_V(Tfinal, x[i]);
    //     double a,b,c;
    //     manufactured_u(Tfinal, x[i], a, b, c);
    //     u1e[i]=a; u2e[i]=b; u3e[i]=c;
    // }
    // auto eV  = error_norms(Vm,   Vex);
    // auto eU1 = error_norms(u1,  u1e);
    // auto eU2 = error_norms(u2,  u2e);

    // RunResult r;
    // r.N = N; r.dx = 1.0/N; r.dt = dt; r.steps = 1; //nsteps;
    // r.errV_inf  = eV.Linf;
    // r.errU1_inf = eU1.Linf;
    // r.errU2_inf = eU2.Linf;

    // Info<< "N   dx   dt   steps      errV_inf"
    //     << "      errU1_inf     errU2_inf\n";

    // Info<< r.N << " "
    //     << r.dx << " "
    //     << r.dt << " "
    //     << r.steps << " "
    //     << r.errV_inf << " "
    //     << r.errU1_inf << " "
    //     << r.errU2_inf
    //     << "\n";

    /*
        // Update the old-time STATES within the ionic model
        refCast<manufacturedFDA>(*ionicModelFDA).updateStatesOld();

        while (pimple.loop())
        {
            // // Before solving the ionic model, reset its state to the old time
            // refCast<manufacturedFDA>(*ionicModelFDA).resetStatesToStatesOld();

            // // Solve the ionic model for the manufactured solutions with 3 gating variable.
            // scalarField& ionicCurrentI = ionicCurrent;
            // scalarField& u1I = u1;
            // scalarField& u2I = u2;
            // scalarField& u3I = u3;

            // // DO NOT reset values
            // // ionicCurrentI = 0.0;
            // // u1I = 0.0;
            // // u2I = 0.0;
            // // u3I = 0.0;
            // ionicModelFDA->calculateCurrent
            // (
            //     runTime.value() - runTime.deltaTValue(),
            //     runTime.deltaTValue(),
            //     Vm.internalField(),
            //     ionicCurrentI,
            //     u1I,
            //     u2I,
            //     u3I
            // );
            // Vm.correctBoundaryConditions();
            // u1.correctBoundaryConditions();
            // u2.correctBoundaryConditions();
            // u3.correctBoundaryConditions();

            // TESTING - enforce exact u1, u2, u3, and ionicCurrent
            // computeManufacturedU
            // (
            //     u1.primitiveFieldRef(),
            //     u2.primitiveFieldRef(),
            //     u3.primitiveFieldRef(),
            //     mesh.C().component(vector::X),
            //     runTime.value()
            // );
            // u1.correctBoundaryConditions();
            // u2.correctBoundaryConditions();
            // u3.correctBoundaryConditions();
            // computeIion
            // (
            //     ionicCurrentI,
            //     u1I,
            //     u2I,
            //     u3I,
            //     Vm,
            //     Cm.value(),
            //     -1.1, // beta
            //     chi.value()
            // );
            // ionicCurrent.correctBoundaryConditions();

            // Construct and solve the voltage equation given a known ionic current
            // and  external stimulus current
            fvScalarMatrix VmEqn
            (
                chi*Cm*fvm::ddt(Vm)
             ==
                fvm::laplacian(conductivity, Vm)
              - chi*ionicCurrent
            );

            VmEqn.solve();

            // Update ionic model using 1st order Euler
            {
                const scalar dt = runTime.deltaTValue();
                const scalarField d
                (
                    u1.primitiveField() + u3.primitiveField() - Vm.primitiveField()
                );
                const scalarField VmMinusU3
                (
                    Vm.primitiveField() - u3.primitiveField()
                );
                const scalarField u2_2(sqr(u2.primitiveField()));
                const scalarField u2_3(u2_2 * u2.primitiveField());

                const scalarField u1new
                (
                    u1.primitiveField()
                  + dt*(sqr(d)*u2.primitiveField() + 0.5*d*u2_2*VmMinusU3)
                );
                const scalarField u2new(u2.primitiveField() + dt*(-d*u2_3));

                // under-relaxation on the ODE update (optional)
                u1.primitiveFieldRef() = u1new;
                u2.primitiveFieldRef() = u2new;
                // u1.primitiveFieldRef() = (1.0-urf)*u1.primitiveField() + urf*u1new;
                // u2.primitiveFieldRef() = (1.0-urf)*u2.primitiveField() + urf*u2new;

                u1.correctBoundaryConditions();
                u2.correctBoundaryConditions();
            }
        }
            */

        // advanceOneStepPicard
        // (
        //     Vm,
        //     u1,
        //     u2,
        //     u3,
        //     Cm,
        //     chi,
        //     -1.1, //beta,
        //     dimensionedScalar
        //     (
        //         "conductivityScalar",
        //         pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
        //         1.1/sqr(constant::mathematical::pi) //electroActivationProperties
	//     ),
        //     runTime.deltaTValue()
        // );

        // runTime.write();

        // computeAndPrintErrors
        // (
        //     Vm, u1, u2, u3, mesh.C().component(vector::X), runTime.value()
        // );
    // }

    // Info<< nl << endl;

    // runTime.printExecutionTime(Info);

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
