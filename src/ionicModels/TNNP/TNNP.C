
/*---------------------------------------------------------------------------* \
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
#include "TNNP_2004.H"
#include "TNNP.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
//#include "Switch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TNNP, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, TNNP, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TNNP::TNNP
(
    const dictionary& dict,
     const label num, 
     const scalar initialDeltaT, 
     const Switch solveVmWithinODESolver
)
:
    ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    STATES_OLD_(num),
    CONSTANTS_(50, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{
    


    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(17, 0.0));
        ALGEBRAIC_.set(i, new scalarField(69, 0.0));
        RATES_.set(i, new scalarField(17, 0.0));

        const dictionary& stimulus = dict.subDict("stimulusParameters");

        TNNPinitConsts
        (
            CONSTANTS_.data(), RATES_[i].data(), STATES_[i].data(), tissue(), stimulus
        );
        {
            // --- Read stimulus parameters from the dictionary ---
            if (!stimulus.found("stim_start")
             || !stimulus.found("stim_period_S1")
             || !stimulus.found("stim_duration")
             || !stimulus.found("stim_amplitude"))
            {
                FatalErrorInFunction
                    << "Missing entries in 'stimulusParameters' sub-dictionary"
                    << abort(FatalError);
            }
    
            CONSTANTS_[5] = readScalar(stimulus.lookup("stim_start"));
            CONSTANTS_[6] = readScalar(stimulus.lookup("stim_period_S1"));
            CONSTANTS_[7] = readScalar(stimulus.lookup("stim_duration"));
            CONSTANTS_[8] = readScalar(stimulus.lookup("stim_amplitude"));
            CONSTANTS_[46]  = readScalar(stimulus.lookup("nstim1"));         
            CONSTANTS_[47] = readScalar(stimulus.lookup("stim_period_S2")); 
            CONSTANTS_[48] = readScalar(stimulus.lookup("nstim2"));         
        }
    }
    Info<< nl << "Initialize ODE constants: " << endl;
    Info<< nl << CONSTANTS_ << endl;
}
    



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TNNP::~TNNP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::TNNP::calculateCurrentSC
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& totalJ
)
{

    const label nIntegrationPoints = STATES_.size();
    label monitorCell = 0;

    const scalar tStart = stepStartTime*1000;
    const scalar tEnd = (stepStartTime + deltaT)*1000;

    forAll(STATES_, integrationPtI)
    {

        // Take a reference to the variables for this integration point
        scalarField& STATESI = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI = RATES_[integrationPtI];


        if (!solveVmWithinODESolver())
        {
            // Vm is evolved by the ODE solver
            STATESI[0] = Vm[integrationPtI]*1000;
        }
        else
        {
            // Set step to deltaT
            scalar& step = ionicModel::step()[integrationPtI];

            step = min(step, deltaT*1000);

            // Update ODE system
            odeSolver().solve(tStart, tEnd, STATESI, step);

            // Calculate the three currents
            ::TNNPcomputeVariables
            (
                tEnd,
                CONSTANTS_.data(),
                RATESI.data(),
                STATESI.data(),
                ALGEBRAICI.data(),
                tissue(),
                solveVmWithinODESolver()
            );
            if (integrationPtI == monitorCell)
            {
                Info<< "integrationPtI = " << integrationPtI
                    << " | t = " << tStart
                    << " → " << tEnd
                    << " | Vm = " << STATESI[0]
                    << " | Iion = " << ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
                    + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
                    + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60]
                    << endl;
            }

            // Extract the total ionic current
            totalJ[integrationPtI] =
                ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
            + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
            + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];
        }
    }
}










void Foam::TNNP::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(69, 0.0);

    // Calculate the rates using the cellML header file
    ::TNNPcomputeRates
    (
        t,
        CONSTANTS_.data(),
        // RATES.data(),
        // STATES.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALGEBRAIC_TMP.data(),
        tissue(), 
        solveVmWithinODESolver()
    );
}


void Foam::TNNP::writeHeader(OFstream& output) const
{
    output << "Time Vm";

    for (int i = 0; i < 17; ++i)
        output << " " << TNNP_STATES_NAMES[i];

    for (int i = 0; i < 69; ++i)
        output << " " << TNNP_ALGEBRAIC_NAMES[i];

    for (int i = 0; i < 17; ++i)
        output << " RATES_" << TNNP_STATES_NAMES[i];

    output << endl;
}


void Foam::TNNP::write(const scalar t, OFstream& output) const

{

    output.precision(8);     
    output.setf(std::ios::scientific);
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


// ************************************************************************* //
