
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
#include "tentusscher_noble_noble_panfilov_2004.H"
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
    const dictionary& dict, const label num, const scalar initialDeltaT
)
:
    ionicModel(dict, num, initialDeltaT),
    STATES_(num),
    CONSTANTS_(47, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{
    // Read the epicardium flag
    //const Switch epicardiumFlag(dict.lookup("epicardium"));

    // Create the integration point lists
    Info<< nl << "Calling initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(17, 0.0));
        ALGEBRAIC_.set(i, new scalarField(69, 0.0));
        RATES_.set(i, new scalarField(17, 0.0));

        // Initialise the constants (repeatedly! it's ok...) and the rates and
        // states
        TNNPinitConsts
        (
            CONSTANTS_.data(), RATES_[i].data(), STATES_[i].data(), tissue()
        );
    }

    if (debug)
    {
        Info<< nl
            << "CONSTANTS = " << CONSTANTS_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TNNP::~TNNP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::TNNP::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& totalJ
)
{
    const label nIntegrationPoints = STATES_.size();

    if (totalJ.size() != nIntegrationPoints)
    {
        FatalErrorInFunction
            << "totalJ.size() != nIntegrationPoints" << abort(FatalError);
    }

    if (Vm.size() != nIntegrationPoints)
    {
        FatalErrorInFunction
            << "Vm.size() != nIntegrationPoints" << abort(FatalError);
    }

    // Update the ODE system for each integration point
    // TODO: this only makes sense if the time-stpe is solved once: otherwise I
    // need to store the old values and only update them for new time-steps
    forAll(STATES_, integrationPtI)
    {
        // Info<< "integrationPtI = " << integrationPtI << endl;

        // Take a reference to the variables for this integration point
        scalarField& STATESI = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI = RATES_[integrationPtI];

        // Update the voltage
        STATESI[0] = Vm[integrationPtI]*1000;

        // Truncate the fields so voltage is ignored
        // scalarField STATES_TRUNCATEDI(STATESI.size() - 1, 0);
        // for (int i = 1; i < STATESI.size(); ++i)
        // {
        //     STATES_TRUNCATEDI[i - 1] = STATESI[i];
        // }

        // Info<< "STATES = " << STATESI << nl
        //     << "ALGEBRAICI = " << ALGEBRAICI << nl
        //     << "RATESI = " << RATESI << endl;
            // << "STATES_TRUNCATEDI = " << STATES_TRUNCATEDI << endl;
            
        // Set t to the initial time and step size
        // Note: ODE solves updates these so we re-create them for each point
        // scalar t = stepStartTime;
        const scalar tStart = stepStartTime*1000;
        const scalar tEnd = (stepStartTime + deltaT)*1000;

         // Set step to deltaT
        scalar& step = ionicModel::step()[integrationPtI];

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
            tissue()
        );

        // Extract the total ionic current
        totalJ[integrationPtI] =
            ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
          + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
          + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];
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
        tissue()
    );
}


// ************************************************************************* //
