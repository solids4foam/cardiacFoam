/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::NashPanfilovActiveTensionModel

Description
    Calculate active muscle tension/stress in terms of the transmembrane
    potential following the approach of MP Nash and AV Panfilov (2004)
    Electromechanical model of excitable tissue to study reentrant cardiac
    arrhythmias, Prog. Biophys. Mol. Biol. 85(2) 501–522.
    
    The active tension, Ta, is given by a 1st-order ordinary differential
    equation

        d(Ta)/dt = epsilonU*(kTa*u - Ta)
    
    where
        epsilonU = 1     if u > 0
        epsilonU = 10    if u < 0
        u = (Vm + 84)/85.7    with Vm in mV
        kTa is constant with dimensions of stress

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef NashPanfilovActiveTensionModel_H
#define NashPanfilovActiveTensionModel_H

#include "ODESystem.H"
#include "dictionary.H"
#include "volFields.H"
#include "ODESolver.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NashPanfilovActiveTensionModel Declaration
\*---------------------------------------------------------------------------*/

class NashPanfilovActiveTensionModel
:
    public ODESystem
{
    // Private Data

        //- Const reference to the dict
        const dictionary& dict_;

        //- Number of integration points, e.g. number of cells
        //  One ODE system will be created for each integration point
        const label nIntegrationPoints_;

        //- Index for monitor point
        //  Gating variables for this point are written to the monitor file
        const label monitorID_;

        //- Monitor file
        OFstream monitorFile_;

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- Internal variables for each integration point: the first entry is
        //  the active tension and the second is the normalised transmembrane
        //  potential
        PtrList<scalarField> internalVariables_;

        //- Model parameter
        const scalar kTa_;

public:

    //- Runtime type information
    TypeName("BuenoOrovio");


    // Constructors

        //- Construct from a dictionary and the number of integration points
        NashPanfilovActiveTensionModel(const dictionary& dict, const label num);


    //- Destructor
    virtual ~NashPanfilovActiveTensionModel()
    {}


    // Member Functions

        //- Solve and ODE system for the start time, the time-step and the
        //  voltage. Then place the active tension in the integration point
        //  result list
        virtual void calculateTension
        (
            const scalar stepStartTime,
            const scalar timeStep,
            const scalarField& Vm,
            scalarField& result
        );

        //- Number of ODE's to solve
        inline virtual label nEqns() const
        {
            return 1;
        }

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const;

        //- Jacobian for solving the coupled system of ODEs
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
