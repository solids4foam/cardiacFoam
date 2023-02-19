/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "NashPanfilovActiveTensionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NashPanfilovActiveTensionModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::NashPanfilovActiveTensionModel::NashPanfilovActiveTensionModel
(
    const dictionary& dict,
    const label num
)
:
    ODESystem(),
    dict_(dict),
    nIntegrationPoints_(num),
    monitorID_(readInt(dict.lookup("monitorID"))),
    monitorFile_("monitorFile"),
    odeSolver_(ODESolver::New(*this, dict_)),
    internalVariables_(num),
    kTa_(readScalar(dict.lookup("kTa")))
{
    // Lookup initial Ta
    const scalar initialTa
    (
        dict.lookupOrDefault<scalar>("initialTa", 0.0143524)
    );

    // Initialise the internal variables
    forAll(internalVariables_, integrationPtI)
    {
        // Two internal variables: Ta and u
        internalVariables_.set(integrationPtI, new scalarField(nEqns() + 1, 0.0));
        internalVariables_[integrationPtI][0] = initialTa;
    }

    // Check monitor is within range
    if (monitorID_ < 0 && monitorID_ > (nIntegrationPoints_ - 1))
    {
        FatalError
            << "monitorID is out of range!" << abort(FatalError);
    }

    // Write header to monitor file
    monitorFile_
        << "# t step Ta u" << endl;
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::NashPanfilovActiveTensionModel::calculateTension
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& activeTension
)
{
    if (activeTension.size() != nIntegrationPoints_)
    {
        FatalErrorIn
        (
            "void Foam::NashPanfilovActiveTensionModel::calculateCurrent(...)"
        )   << "activeTension.size() != nIntegrationPoints" << abort(FatalError);
    }

    if (Vm.size() != nIntegrationPoints_)
    {
        FatalErrorIn
        (
            "void Foam::NashPanfilovActiveTensionModel::calculateCurrent(...)"
        )   << "Vm.size() != nIntegrationPoints" << abort(FatalError);
    }

    // Update the ODE system for each integration point
    // TODO: this only makes sense if the time-step is solved once: otherwise I
    // need to store the old values and only update them for new time-steps
    // One solution would be to use volScalarFields; an added benefit is
    // visualisation, but a downside is that we cannot easily use pointFields or
    // surfaceFields
    forAll(internalVariables_, integrationPtI)
    {
        // Take a reference to the internal variables for this integration point
        scalarField& yStart = internalVariables_[integrationPtI];

        // Set t to the initial time and step size
        // Note: ODE solves updates these so we re-create them for each point
        const scalar tStart = stepStartTime;
        const scalar tEnd = stepStartTime + deltaT;
        scalar step = deltaT;

        // Take references to improve readability
        const scalar& Ta = yStart[0];
        scalar& u = yStart[1];

        // Calculate the normalised voltage
        u = (1000*Vm[integrationPtI] + 84)/85.7;

        // Update ODE system
        // Note: this updates the value of t to the end time
        odeSolver_->solve(tStart, tEnd, yStart, step);

        // Retrieve the active tension
        activeTension[integrationPtI] = Ta;

        // Write to the monitor file
        if (integrationPtI == monitorID_)
        {
            monitorFile_
                << tEnd << " "
                << step << " "
                << Ta << " "
                << u << endl;
        }
    }
}


void Foam::NashPanfilovActiveTensionModel::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Note: it is assumed that u_ holds the voltage for the current integration point

    // Take references to improve readability
    const scalar& Ta = y[0];
    const scalar& u = y[1];
    scalar& dTadt = dydt[0];

    // EpsilonU parameter
    scalar epsilonU = 0;
    if (u > 0)
    {
        epsilonU = 1.0;
    }
    else
    {
        epsilonU = 10.0;
    }

    // Calculate derivative
    dTadt = epsilonU*(kTa_*u - Ta);    
}


void Foam::NashPanfilovActiveTensionModel::jacobian
(
    const scalar t,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("void Foam::NashPanfilovActiveTensionModel::jacobian(...)");
}


// ************************************************************************* //
