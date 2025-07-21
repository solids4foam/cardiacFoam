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
#include "Courtemanche_1998.H"
#include "Courtemanche.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"

//Only needs strings for the header writing 
//#include <string>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Courtemanche, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, Courtemanche, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::Courtemanche::Courtemanche
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

    // Create the integration point lists
    Info<< nl << "Calling Courtemanche initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(NUM_STATES, 0.0));
        ALGEBRAIC_.set(i, new scalarField(NUM_CONSTANTS, 0.0));
        RATES_.set(i, new scalarField(NUM_STATES, 0.0));

        // Initialise the constants (repeatedly! it's ok...) and the rates and
        // states
        CourtemancheinitConsts
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

Foam::Courtemanche::~Courtemanche()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::Courtemanche::calculateCurrent
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
            //STATESI[0] = Vm[integrationPtI]*1000;
            STATESI[membrane_V] = Vm[integrationPtI]*1000; // 
            
            //???? Ask Philip. No need to define in the temporary saver,
            // it can be state[0] still I think and it maps membrane_v


            const scalar tStart = stepStartTime*1000;
            const scalar tEnd = (stepStartTime + deltaT)*1000;

            // Set step to deltaT
            scalar& step = ionicModel::step()[integrationPtI];

            // Update ODE system
            odeSolver().solve(tStart, tEnd, STATESI, step);

            // Calculate the three currents
            ::CourtemanchecomputeVariables
            (
                tEnd,
                CONSTANTS_.data(),
                RATESI.data(),
                STATESI.data(),
                ALGEBRAICI.data(),
                tissue()
            );

        // Extract the total ionic current
        totalJ[integrationPtI] =ALGEBRAICI[AV_INa] + ALGEBRAICI[AV_IK1] + ALGEBRAICI[AV_Ito] + 
        ALGEBRAICI[AV_IKur] + ALGEBRAICI[AV_IKr] + ALGEBRAICI[AV_IKs] + ALGEBRAICI[AV_ICaL] + 
        ALGEBRAICI[AV_IpCa] + ALGEBRAICI[AV_INaK] + ALGEBRAICI[AV_INaCa] + ALGEBRAICI[AV_IbNa] + 
        ALGEBRAICI[AV_IbCa];
    


    }
}

void Foam::Courtemanche::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // VOI is the curren time
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    // Calculate the rates using the cellML header file
    ::CourtemanchecomputeVariables
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



// void Foam::Courtemanche::jacobian
// (
//     const scalar t,
//     const scalarField& y,
//     scalarField& dfdt,
//     scalarSquareMatrix& dfdy
// ) const
// {
//     notImplemented("void jacobian(...)");
// }

// ************************************************************************* //


/*---------------------------------------------------------------------------*\
This two functions are related to the single cell solver. I still need to undrstand the best way to implement it 
but I would say, since it has a slighly diferent solver.C file, it will have a different src?
I dontkn know because the models are the same, maybe just the ionic model files differ. 


void Foam::ionicModelCellML::writeHeader(OFstream& output) const
{
    output << "time";

    for (int i = 0; i < NUM_STATES; ++i)
    {
        output << " " << STATES_NAMES[i];
    }
    for (int i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        output << " " << ALGEBRAIC_NAMES[i];
    }
    for (int i = 0; i < NUM_STATES; ++i)
    {
        output << "RATES_ " << STATES_NAMES[i];
    }

    output 
        << endl;
}
void Foam::ionicModelCellML::write(const scalar t, OFstream& output) const
{
    // Write the results
    Info<< "Writing to " << output.name() << endl;
    output
        << t;
    forAll(STATES_, i)
    {
        output
            << " " << STATES_[i];
    }
    forAll(ALGEBRAIC_, i)
    {
        output
            << " " << ALGEBRAIC_[i];
    }
    forAll(RATES_, i)
    {
        output
            << " " << RATES_[i];
    }
    output
        << endl;
}

\*---------------------------------------------------------------------------*/

