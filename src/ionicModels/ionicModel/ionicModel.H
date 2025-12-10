/*---------------------------------------------------------------------------*\
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

Class
    Foam::ionicModel

Description
    Abstract base class for ionic models.

SourceFiles
    ionicModel.C

Authors
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ionicModel_H
#define ionicModel_H

#include "runTimeSelectionTables.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "word.H"
#include "Switch.H"
#include "OFstream.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ionicModel Declaration
\*---------------------------------------------------------------------------*/

class ionicModel
:
    public ODESystem
{
    // Private Data

    // Private Data

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- Copy of the model dictionary
        dictionary dict_;

        //- Step size for ODE solver
        mutable scalarList step_;

        //- Flag for tissue type
        mutable label tissue_;

        //- Flag for enable the update of the transmembrane voltage within the
        //  the ionic model ODE solver. For example, if calling for a single
        //  cell, then the governing voltage PDE becomes an ODE so you would
        //  enable this flag
        mutable Switch solveVmWithinODESolver_;


    // Private Member Functions

        //- No copy construct
        ionicModel(const ionicModel&) = delete;

        //- No copy assignment
        void operator=(const ionicModel&) = delete;

protected:
        virtual void setTissueFromDict();
        wordList variableExport_;

        // // Debug settings
        // bool debugMode_ = false;
        // label monitorCell_ = -1;
        wordList debugVarNames_;

        /*NOTE: currently protected, each means they are only visible by the derived class.
        THe higher order solver will use them so probably will be set as public.*/

        // Copies STATES_ → states external buffer
        inline void copyInternalToExternal
        (
            const PtrList<scalarField>& internal,
            Field<Field<scalar>>& external,
            const label nStates
        ) const
        {
            forAll(internal, i)
            {
                if (external[i].size() >= nStates)
                    external[i] = internal[i];
            }
        }

        // Copies current STATES_ → STATES_OLD_
        inline void saveStateSnapshot
        (
            const PtrList<scalarField>& internal,
            PtrList<scalarField>& old
        ) const
        {
            forAll(internal, i)
                old[i] = internal[i];
        }

        // Restores snapshot STATES_OLD_ → STATES_
        inline void restoreStateSnapshot
        (
            PtrList<scalarField>& internal,
            const PtrList<scalarField>& old
        ) const
        {
            forAll(internal, i)
                internal[i] = old[i];
        }


public:

    //- Runtime type information
    TypeName("ionicModel");


    // Declare run-time constructor selection table

        //- Ionic model run-time selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ionicModel,
            dictionary,
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver
            ),
            (dict, nIntegrationPoints, initialDeltaT, solveVmWithinODESolver)
        );

    // Selectors

        //- Create a ionic model
        static autoPtr<ionicModel> New
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );

    // Constructors

        //- Construct given a dictionary, the number of integration points and a
        //  initial time step for the ODE solver
        ionicModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );



    //- Destructor
    virtual ~ionicModel();


//-----------------------Member Functions ---------------------//

//-----------------------Solver Functions ---------------------//
        //- ODE solver
        ODESolver& odeSolver()
        {
            if (!odeSolver_)
                {odeSolver_.reset(ODESolver::New(*this, dict_));}
            return odeSolver_();
        }

        //- Step size for ODE solver
        scalarList& step()
            {return step_;}

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        )= 0;

        //- Number of ODE's to solve
        virtual label nEqns() const = 0;

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation - right side of the equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const = 0;


        //- Solve and ODE system for the start time, the time-step and the
        //  voltage. Then place the result in the integration point result list
        virtual void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        )= 0;

        //Initialization for a volumetric field
        // - currently for manufatured solutions
        virtual void initializeFields
        (
            volScalarField& Vm,
            volScalarField& u1,
            volScalarField& u2,
            volScalarField& u3,
            const volVectorField& C
        )
        {
        }

        //- Jacobian for solving the coupled system of ODEs
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const
        {
            notImplemented
            (
                "Please define this function in the derived class!"
            );
        };


//-----------------------Read and get for derived class ---------------------//

        //- Getter for tissue flag
        label tissue() const
        {
            return tissue_;
        }

        //- Getter for 1 cell inside the CellML ODE
        const Switch& solveVmWithinODESolver() const
        {
            return solveVmWithinODESolver_;
        }

        //Tissue and dimensions flag getter for ionicModel and manufactured solution
            //tissue
        virtual List<word> supportedTissueTypes() const
        {
            return {"epicardialCells", "mCells", "endocardialCells", "myocyte"};
        }

        virtual word tissueName() const
        {
        const label t = tissue_;

        return (t == 1) ? "epicardialCells" :
            (t == 2) ? "mCells" :
            (t == 3) ? "endocardialCells" :
            (t == 4) ? "myocyte" :
                        "unknown";
        }
            //dimension in manufactured
        virtual bool hasManufacturedSolution() const
        {
        return false;   // Default: no manufactured solution
        }

        virtual List<word> supportedDimensions() const
        {
            return {"1D", "2D", "3D"};
        }


//-------------------------------Print Export and write ---------------------//
        // Write and export from internal variables functions
        // singleCell
        virtual void writeHeader(OFstream& os) const
        {}

        virtual void write(const scalar t, OFstream& os) const
        {}

        // 3D vol fields(defined for states, algebraic and rates.)
        virtual void exportStates
        (
            const Field<Field<scalar>>&,
            PtrList<volScalarField>&
        )
        {}

        // Printing function
        virtual void debugPrintFields
        (
            label,
            scalar,
            scalar,
            scalar
        ) const
        {
        }

        // Return names of variables to export
        virtual wordList exportedFieldNames() const
        {
            return variableExport_.size()
                ? variableExport_
                : wordList();
        }
        virtual wordList debugPrintedNames() const
        {
            return debugVarNames_.size()
                ? debugVarNames_
                : wordList();
        }

    };

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

