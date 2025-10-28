/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include <math.h>
#include "manufacturedFDA.H"
#include "manufacturedFDA_2014.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModelFDA.H"
#include "fieldInit.H"
//Only needs strings for the header writing
//#include <string>

#include "manufacturedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manufacturedFDA, 0);
    addToRunTimeSelectionTable
    (
        ionicModelFDA, manufacturedFDA, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manufacturedFDA::manufacturedFDA
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ionicModelFDA(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    STATES_OLD_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{

    //see if I need to add flog in function as well.
    Info<< nl << "Calling FDA test Constants" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(NUM_STATES, 0.0));
        STATES_OLD_.set(i, new scalarField(NUM_STATES, 0.0));
        ALGEBRAIC_.set(i, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(i, new scalarField(NUM_STATES, 0.0));



        // Initialise the constants (repeatedly! it's ok...) and the rates and
        // states
        manufacturedFDAinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue()
        );

        // Initialise old time STATES
        STATES_OLD_[i] = STATES_[i];
    }

    Info<< nl
        << "CONSTANTS = " << CONSTANTS_ << endl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manufacturedFDA::~manufacturedFDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::manufacturedFDA::supportedTissues() const
{
    return {"1D", "2D", "3D"};
}



 // --- Initialize OpenFOAM fields and ODE states ---
void Foam::manufacturedFDA::initializeFields(
    volScalarField& Vm,
    volScalarField& u1m,
    volScalarField& u2m,
    volScalarField& u3m,
    const volVectorField& C
)
{
    // Initialize OpenFOAM fields
    forAll(C, celli)
    {

        const vector& coord = C[celli];
        Info << "C size: " << C.size() << ", Vm size: " << Vm.size() << nl;

        double x = coord.x();
        double y = coord.y();
        double z = coord.z();

        double VmVal, uu1, uu2, uu3;
        u_init(x, y, z, VmVal, uu1, uu2, uu3, tissue());

        Vm[celli]  = VmVal;
        u1m[celli] = uu1;
        u2m[celli] = uu2;
        u3m[celli] = uu3;

        Info << "Cell " << celli
             << " x=" << x << " y=" << y << " z=" << z
             << " Vm=" << Vm[celli]
             << " u1=" << u1m[celli]
             << " u2=" << u2m[celli]
             << " u3=" << u3m[celli]
             << endl;
    }

    // Initialize STATES_ array for ODE solver
    forAll(STATES_, i)
    {
        scalarField& STATESI = STATES_[i];
        STATESI[V]  = Vm[i];
        STATESI[u1] = u1m[i];
        STATESI[u2] = u2m[i];
        STATESI[u3] = u3m[i];

        // Copy to old-time STATES
        STATES_OLD_[i] = STATESI;
    }
    //correct boundary conditions
    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    Info << "Boundary conditions corrected for Vm, u1m, u2m, and u3m." << endl;
}




void Foam::manufacturedFDA::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    scalarField& u1m,
    scalarField& u2m,
    scalarField& u3m
)
{
    const label nIntegrationPoints = STATES_.size();

    if (Im.size() != nIntegrationPoints)
    {
        FatalErrorInFunction
            << "I.size() != nIntegrationPoints" << abort(FatalError);
    }

    if (Vm.size() != nIntegrationPoints)
    {
        FatalErrorInFunction
            << "Vm.size() != nIntegrationPoints" << abort(FatalError);
    }

    // Reset states to their old time values
    // This allows the use of PIMPLE-type outer iterations within the solver,
    // i.e. this ODE solver can then potentially be called multiple times in the
    // same time step
    Info<< "STATES_ = " << STATES_ << endl;

    STATES_ = STATES_OLD_;

    Info<< "STATES_ = " << STATES_ << endl;

    // Debug - calculate true solutions
    // I don't have access to the mesh so I will hard-code the coordinates
    // scalarField exactU1(STATES_.size());
    // scalarField exactU2(STATES_.size());
    // scalarField exactU3(STATES_.size());
    // scalarField x(STATES_.size());
    // const scalar dx = 1.0/scalar(STATES_.size());
    // forAll(x, i)
    // {
    //     x[i] = (scalar(i) + 0.5)*dx;
    // }
    // computeManufacturedU(exactU1, exactU2, exactU3, x, stepStartTime + deltaT);
    // scalarField exactIion(STATES_.size());;
    // computeIion
    // (
    //     exactIion, exactU1, exactU2, exactU3, Vm,
    //     CONSTANTS_[Cm], CONSTANTS_[Beta], CONSTANTS_[Chi]
    // );



    label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        // Take a reference to the variables for this integration point
        scalarField& STATESI = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI = RATES_[integrationPtI];

        const scalar tStart = stepStartTime;
        const scalar tEnd = (stepStartTime + deltaT);

        // Set step to deltaT
        scalar& step = ionicModelFDA::step()[integrationPtI];

        // Update the voltage used by the ODE solver - THIS WAS MISSING
        STATESI[V] = Vm[integrationPtI];

        //clamp the step, or just leave it to the solver
        step = min(step, deltaT);

        if (integrationPtI == monitorCell)
        {
            Info<< "integrationPtI = " << integrationPtI
                << " | t = " << tStart
                << " → " << tEnd
                << " | step = " << step
                << " | Vm = " << STATESI[V]
                << " | u1 = " << STATESI[u1]
                << " | u2 = " << STATESI[u2]
                << " | u3 = " << STATESI[u3]
                << endl;
        }

        // Update ODE system
        odeSolver().solve(tStart, tEnd, STATESI, step);

        // Calculate the three currents
        ::manufacturedFDAcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        Im[integrationPtI] = ALGEBRAICI[Iion];
        u1m[integrationPtI]= STATESI[u1];
        u2m[integrationPtI]= STATESI[u2];
        u3m[integrationPtI]= STATESI[u3];
    }
}


void Foam::manufacturedFDA::writeHeader(OFstream& output) const
{

    output << "time Vm";

    for (int i = 0; i < NUM_STATES; ++i)
        output << " " << manufacturedFDASTATES_NAMES[i];

    for (int i = 0; i < NUM_ALGEBRAIC; ++i)
        output << " " << manufacturedFDAALGEBRAIC_NAMES[i];

    for (int i = 0; i < NUM_STATES; ++i)
        output << " RATES_" << manufacturedFDASTATES_NAMES[i];

    output << endl;
}


void Foam::manufacturedFDA::write(const scalar t, OFstream& output) const

{

    scalar Vm = STATES_[0][0];

    output
        << t << " " << Vm;

    // States
    forAll(STATES_[0], j)
    {
        output << " " << STATES_[0][j];
    }

    // Algebraic variables
    forAll(ALGEBRAIC_[0], j)
    {
        output << " " << ALGEBRAIC_[0][j];
    }

    // Rates
    forAll(RATES_[0], j)
    {
        output << " " << RATES_[0][j];
    }


output << endl;
}


void Foam::manufacturedFDA::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(1, 0.0);

    compute_f_RHS_point
    (
        dydt[1],
        dydt[2],
        dydt[3],
        y[1],
        y[2],
        y[3],
        y[0]
    );

    // // Calculate the rates using the cellML header file
    // ::manufacturedFDAcomputeVariables
    // (
    //     t,
    //     CONSTANTS_.data(),
    //     // RATES.data(),
    //     // STATES.data(),
    //     dydt.data(),
    //     const_cast<scalarField&>(y).data(),
    //     ALGEBRAIC_TMP.data(),
    //     tissue(),
    //     solveVmWithinODESolver()
    // );
}
