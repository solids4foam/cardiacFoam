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
#include "tentusscher_noble_noble_panfilov_2004.H"
#include "ionicModelCellML.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ionicModelCellML, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModelCellML::ionicModelCellML
(
    const dictionary& dict, const label num
)
:
    ODESystem(),
    odeSolver_(ODESolver::New(*this, dict)),
    STATES_(num),
    CONSTANTS_(46, 0.0),
    ALGEBRAIC_(num),
    RATES_(num),
    step_(num, readScalar(dict.lookup("initialODEStep")))
{
    // Create the integration point lists
    Info<< nl << "Calling initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(17, 0.0));
        ALGEBRAIC_.set(i, new scalarField(69, 0.0));
        RATES_.set(i, new scalarField(17, 0.0));

        // Initialise the constants (repeatedly! it's ok...) and the rates and
        // states
        initConsts(CONSTANTS_.data(), RATES_[i].data(), STATES_[i].data());
    }

    if (debug)
    {
        Info<< nl
            << "CONSTANTS = " << CONSTANTS_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModelCellML::~ionicModelCellML()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::ionicModelCellML::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& totalJ
)
{
    // Info<< "STATES_[0] " << STATES_[1] << nl
    //     << "STATES_[1] " << STATES_[1] << nl
    //     << "STATES_[2] " << STATES_[1] << endl;

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
        scalar& step = step_[integrationPtI];

        // Info<< "tStart = " << tStart << nl
        //     << "tEnd = " << tEnd << nl
        //     << "step = " << step << endl;

        // Update ODE system
        // Info<< "odeSolver_->solve" << endl;
        // odeSolver_->solve(tStart, tEnd, STATES_TRUNCATEDI, step);
        odeSolver_->solve(tStart, tEnd, STATESI, step);
        // Info<< "odeSolver_->solve: done" << endl;

        // Calculate the three currents
        ::computeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data()
        );

        // Extract the total ionic current
        totalJ[integrationPtI] =
            ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
          + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
          + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];

    // Info<< "STATES_[0] " << STATES_[1] << nl
    //     << "STATES_[1] " << STATES_[1] << nl
    //     << "STATES_[2] " << STATES_[1] << endl;
        // if (integrationPtI == monitorID_)
        // {
        //     monitorFile_
        //         << tEnd << " "
        //         << step << " "
        //         << u << " "
        //         << yStart[0] << " "
        //         << yStart[1] << " "
        //         << yStart[2] << " "
        //         << Jfi << " "
        //         << Jso << " "
        //         << Jsi << endl;
        // }
    }
}


// void Foam::ionicModelCellML::solve
// (
//     const scalar tOld, const scalar deltaT
// ) const
// {
//     // Define the start time, end time and initial time step
//     const scalar tStart = tOld;
//     const scalar tEnd = tOld + deltaT;

//     // Take a reference to the initial state
//     scalarField& yStart = STATES_;

//     scalar step = min(step_, deltaT);

//     // Solve the ODE system
//     Info<< "step = " << step_ << endl;
//     odeSolver_->solve(tStart, tEnd, yStart, step);

//     // Update step
//     step_ = step;
// }


// void Foam::ionicModelCellML::computeVariables(const scalar t) const
// {
//     ::computeVariables
//     (
//         t,
//         CONSTANTS_.data(),
//         RATES_.data(),
//         STATES_.data(),
//         ALGEBRAIC_.data()
//     );
// }


// void Foam::ionicModelCellML::write(const scalar t, OFstream& output) const
// {
//     // Write the results
//     Info<< "Writing to " << output.name() << endl;
//     output
//         << t;
//     forAll(STATES_, i)
//     {
//         output
//             << " " << STATES_[i];
//     }
//     forAll(ALGEBRAIC_, i)
//     {
//         output
//             << " " << ALGEBRAIC_[i];
//     }
//     forAll(RATES_, i)
//     {
//         output
//             << " " << RATES_[i];
//     }
//     output
//         << endl;
// }


void Foam::ionicModelCellML::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // InfoInFunction
    //     << "start" << endl;

    // scalarField ALGEBRAIC_TMP(ALGEBRAIC_[0].size(), 0.0);
    scalarField ALGEBRAIC_TMP(69, 0.0);

    // // Create temporary lists that are the correct size
    // scalarField STATES(y.size() + 1, 0);
    // scalarField RATES(y.size() + 1, 0);
    // for (int i = 1; i < STATES.size(); ++i)
    // {
    //     STATES[i] = y[i - 1];
    //     RATES[i] = dydt[i - 1];
    // }

    // Calculate the rates using the cellML header file
    // Info<< "computeRates start" << endl;
    ::computeRates
    (
        t,
        CONSTANTS_.data(),
        // RATES.data(),
        // STATES.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALGEBRAIC_TMP.data()
    );
    // Info<< "computeRates end" << endl;

    // // Copy RATES into dydt
    // for (int i = 1; i < STATES.size(); ++i)
    // {
    //     dydt[i - 1] = RATES[i];
    // }    

    // InfoInFunction
    //     << "end" << endl;
 }


void Foam::ionicModelCellML::jacobian
(
    const scalar t,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("void jacobian(...)");
}

// ************************************************************************* //
