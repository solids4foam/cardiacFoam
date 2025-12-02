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
#include "tmanufacturedFDA.H"
#include "tmanufacturedFDA_2014.H"
#include "tmanufacturedFDA_2014Names.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
//Only needs strings for the header writing
//#include <string>

#include "tmanufacturedFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tmanufacturedFDA, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, tmanufacturedFDA, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tmanufacturedFDA::tmanufacturedFDA
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
        tmanufacturedFDAinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue()
        );
        STATES_OLD_[i] = STATES_[i];
    }

    Info<< nl
        << "CONSTANTS = " << CONSTANTS_ << endl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tmanufacturedFDA::~tmanufacturedFDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::tmanufacturedFDA::supportedTissueTypes() const

{
    return {"1D", "2D", "3D"};
}


void Foam::tmanufacturedFDA::initializeFields
(
    volScalarField& Vm,
    volScalarField& u1m,
    volScalarField& u2m,
    volScalarField& u3m,
    const volVectorField& C
)
{
    scalarField X = C.component(vector::X);
    scalarField Y = C.component(vector::Y);
    scalarField Z = C.component(vector::Z);
    const scalar t = 0.0;

    computeManufacturedV(Vm, X, Y, Z, t, tissue());
    computeManufacturedU(u1m, u2m, u3m, X, Y, Z, t, tissue());

    forAll(STATES_, i)
    {
        scalarField& S = STATES_[i];
        S[V]  = Vm[i];
        S[u1] = u1m[i];
        S[u2] = u2m[i];
        S[u3] = u3m[i];

        STATES_OLD_[i] = S;
    }

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    Info<< "Boundary conditions corrected." << endl;
}

void Foam::tmanufacturedFDA::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    STATES_ = STATES_OLD_;
    const scalar tStart = stepStartTime;
    const scalar tEnd   = stepStartTime + deltaT;
    const label monitorCell = 0;

    forAll(STATES_, i)
    {
        scalarField& S = STATES_[i];
        scalarField& A = ALGEBRAIC_[i];
        scalarField& R = RATES_[i];

        S[V] = Vm[i];

        ::tmanufacturedFDAcomputeVariables
        (
            tStart,
            CONSTANTS_.data(),
            R.data(),
            S.data(),
            A.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        if (i == monitorCell)
        {
            Info<< "integrationPtI = " << i
                << " | t = " << tStart
                << " → " << tEnd
                << " | Vm = " << S[V]
                << " | u1 = " << S[u1]
                << " | u2 = " << S[u2]
                << " | Iion = " << A[Iion]
                << endl;
        }

        Im[i] = A[Iion];
    }
}

// void Foam::tmanufacturedFDA::calculateGating
// (
//     const scalar stepStartTime,
//     const scalar deltaT,
//     const scalarField& Vm,
//     scalarField& Im,
//     Field<Field<scalar>>& states
// )
// {

// }

void Foam::tmanufacturedFDA::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    STATES_ = STATES_OLD_;
    const scalar tStart = stepStartTime;
    const scalar tEnd   = tStart + deltaT;
    const label monitorCell = 0;

    forAll(STATES_, i)
    {
        scalarField& S = STATES_[i];
        scalarField& A = ALGEBRAIC_[i];
        scalarField& R = RATES_[i];

        scalar& h = ionicModel::step()[i];
        S[V] = Vm[i];
        h = min(h, deltaT);

        if (i == monitorCell)
        {
            Info<< "solveODE: i=" << i
                << " | t = " << tStart << " → " << tEnd
                << " | step = " << h
                << " | Vm = " << S[V]
                << " | u1 = " << S[u1]
                << " | u2 = " << S[u2]
                << " | Iion = " << A[Iion]
                << endl;
        }

        odeSolver().solve(tStart, tEnd, S, h);

        ::tmanufacturedFDAcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            R.data(),
            S.data(),
            A.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        Im[i] = A[Iion];

        // u1m[i] = S[u1];
        // u2m[i] = S[u2];
        // u3m[i] = S[u3];
    }
}


// ------------------------------------------------------------------------- //
//  Manufactured export
// ------------------------------------------------------------------------- //

void Foam::tmanufacturedFDA::exportManufacturedStates
(
    volScalarField& VmMS,
    volScalarField& u1MS,
    volScalarField& u2MS,
    volScalarField& u3MS
) const
{
    forAll(STATES_, i)
    {
        const scalarField& S = STATES_[i];
        VmMS[i] = S[V];
        u1MS[i] = S[u1];
        u2MS[i] = S[u2];
        u3MS[i] = S[u3];
    }

    VmMS.correctBoundaryConditions();
    u1MS.correctBoundaryConditions();
    u2MS.correctBoundaryConditions();
    u3MS.correctBoundaryConditions();
}
Foam::wordList Foam::tmanufacturedFDA::exportedFieldNames() const
{
    // If user defined a list in "exportedVariables", use it
    if (variableExport_.size())
    {
        return variableExport_;
    }

    // Default: export all manufactured states + algebraic Iion
    Foam::wordList names;

    for (Foam::label i = 0; i < NUM_STATES; ++i)
    {
        names.append(tmanufacturedFDASTATES_NAMES[i]);
    }

    for (Foam::label j = 0; j < NUM_ALGEBRAIC; ++j)
    {
        names.append(tmanufacturedFDAALGEBRAIC_NAMES[j]);
    }

    return names;
}
void Foam::tmanufacturedFDA::exportStates
(
    const Field<Field<scalar>>&,
    List<volScalarField*>& outFields
)
{
    const Foam::wordList names = exportedFieldNames();

    if (outFields.size() != names.size())
    {
        WarningInFunction
            << "exportStates: expected " << names.size()
            << " fields, got " << outFields.size() << endl;
        return;
    }

    // Loop over all cells
    forAll(STATES_, cellI)
    {
        const scalarField& S = STATES_[cellI];      // states
        const scalarField& A = ALGEBRAIC_[cellI];   // algebraics

        // For each requested export variable
        for (Foam::label k = 0; k < names.size(); ++k)
        {
            Foam::volScalarField& fld = *outFields[k];
            const Foam::word& name = names[k];

            bool matched = false;

            // ----- MATCH STATE VARIABLES (V, u1, u2, u3) -----
            for (Foam::label s = 0; s < NUM_STATES; ++s)
            {
                if (name == tmanufacturedFDASTATES_NAMES[s])
                {
                    fld[cellI] = S[s];
                    matched = true;
                    break;
                }
            }

            // ----- MATCH ALGEBRAIC VARIABLES (Iion) -----
            if (!matched)
            {
                for (Foam::label a = 0; a < NUM_ALGEBRAIC; ++a)
                {
                    if (name == tmanufacturedFDAALGEBRAIC_NAMES[a])
                    {
                        fld[cellI] = A[a];
                        matched = true;
                        break;
                    }
                }
            }

            if (!matched)
            {
                WarningInFunction
                    << "Unknown variable '" << name
                    << "' in exportedVariables for tmanufacturedFDA" << endl;
            }
        }
    }

    // Correct boundaries for each exported field
    for (Foam::label k = 0; k < names.size(); ++k)
    {
        outFields[k]->correctBoundaryConditions();
    }
}


// ------------------------------------------------------------------------- //
//  IO
// ------------------------------------------------------------------------- //

void Foam::tmanufacturedFDA::writeHeader(OFstream& os) const
{
    os << "time Vm";

    for (int i=0;i<NUM_STATES;++i)
        os << " " << tmanufacturedFDASTATES_NAMES[i];

    for (int i=0;i<NUM_ALGEBRAIC;++i)
        os << " " << tmanufacturedFDAALGEBRAIC_NAMES[i];

    for (int i=0;i<NUM_STATES;++i)
        os << " RATES_" << tmanufacturedFDASTATES_NAMES[i];

    os << endl;
}

void Foam::tmanufacturedFDA::write(const scalar t, OFstream& os) const
{
    os << t << " " << STATES_[0][0];

    forAll(STATES_[0], j)
        os << " " << STATES_[0][j];

    forAll(ALGEBRAIC_[0], j)
        os << " " << ALGEBRAIC_[0][j];

    forAll(RATES_[0], j)
        os << " " << RATES_[0][j];

    os << endl;
}


void Foam::tmanufacturedFDA::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALG(NUM_ALGEBRAIC, 0.0);

    ::tmanufacturedFDAcomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALG.data(),
        tissue(),
        solveVmWithinODESolver()
    );
}
// ************************************************************************* //