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
#include "fvMesh.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"

#include "electromechanicsFeedbackProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ionicModel Declaration
\*---------------------------------------------------------------------------*/

class ionicModel
:
    public ODESystem,
    public CouplingSignalProvider
{
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

        //- External stimulus protocol (separate from model constants)
        StimulusProtocol stimulus_;


    // Private Member Functions

        //- No copy construct
        ionicModel(const ionicModel&) = delete;

        //- No copy assignment
        void operator=(const ionicModel&) = delete;

    protected:
        virtual void setTissueFromDict();
        void setTissue(label t) { tissue_ = t; }
        void setStimulusProtocolFromDict(const dictionary& dict)
        {
            stimulus_ = stimulusIO::loadStimulusProtocol(dict);
        }
        wordList variableExport_;
        wordList debugVarNames_;
        mutable ionicModelIO::FullPlanCache fullWritePlanCache_;
        mutable ionicModelIO::SelectedMapCache writeSelectedPlanCache_;
        mutable ionicModelIO::SelectedMapCache exportSelectedPlanCache_;
        mutable ionicModelIO::SelectedMapCache exportTransferSelectedPlanCache_;
        mutable ionicModelIO::SelectedMapCache importTransferSelectedPlanCache_;
        mutable ionicModelIO::SelectedMapCache debugSelectedPlanCache_;
        mutable ionicModelIO::ExportedNamesCache exportedNamesCache_;
        mutable ionicModelIO::ExportedNamesCache debugNamesCache_;

        // Optional metadata hooks for generic I/O in the base class.
        // Derived models can override these and then rely on base
        // write/export/debug implementations.
        virtual const PtrList<scalarField>* ioStatesPtr() const
        {
            return nullptr;
        }

        virtual const PtrList<scalarField>* ioAlgebraicPtr() const
        {
            return nullptr;
        }

        virtual const PtrList<scalarField>* ioRatesPtr() const
        {
            return nullptr;
        }

        virtual const scalarField* ioConstantsPtr() const
        {
            return nullptr;
        }

        virtual label ioNumConstants() const
        {
            return 0;
        }

        virtual const char* const* ioConstantNames() const
        {
            return nullptr;
        }

        virtual const char* const* ioStateNames() const
        {
            return nullptr;
        }

        virtual label ioNumStates() const
        {
            return 0;
        }

        virtual const char* const* ioAlgebraicNames() const
        {
            return nullptr;
        }

        virtual label ioNumAlgebraic() const
        {
            return 0;
        }

        virtual ionicModelIO::VmTransform ioVmTransform() const
        {
            return nullptr;
        }

        bool hasIOMetadata() const
        {
            return ioStateNames() && ioNumStates() > 0;
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

    /// --- CouplingSignalProvider interface ---
    bool hasSignal(CouplingSignal s) const override;
    scalar signal(label i, CouplingSignal s) const override;
    

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
            scalarField& Im
        ) = 0;

        //- Number of ODE's to solve
        virtual label nEqns() const override = 0;

        //- Derivatives of each ODE with respect to the primary unknown in that
        //  equation - right side of the equation
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const override = 0;


        //- Jacobian for solving the coupled system of ODEs
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const override
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

        const StimulusProtocol& stimulusProtocol() const
        {
            return stimulus_;
        }

        //- Return available state variable names for this model
        wordList stateVariableNames() const
        {
            if (!hasIOMetadata())
            {
                return wordList();
            }

            const label nStates = ioNumStates();
            if (nStates <= 0)
            {
                return wordList();
            }

            const char* const* names = ioStateNames();
            if (!names)
            {
                return wordList();
            }

            wordList result(nStates);
            for (label i = 0; i < nStates; ++i)
            {
                result[i] = word(names[i]);
            }
            return result;
        }

        //- Return constants values for this model
        scalarField constantVariableValues() const
        {
            const auto* constantsPtr = ioConstantsPtr();
            if (!constantsPtr)
            {
                return scalarField();
            }
            return *constantsPtr;
        }

        //- Return available constant variable names for this model
        wordList constantVariableNames() const
        {
            const label nConstants = ioNumConstants();
            if (nConstants <= 0)
            {
                return wordList();
            }

            const char* const* names = ioConstantNames();
            if (!names)
            {
                return wordList();
            }

            wordList result(nConstants);
            for (label i = 0; i < nConstants; ++i)
            {
                result[i] = word(names[i]);
            }
            return result;
        }

        //- Return state values at an integration point
        scalarField stateVariableValues(const label integrationPtI = 0) const
        {
            const auto* statesPtr = ioStatesPtr();
            if (!statesPtr || statesPtr->empty())
            {
                return scalarField();
            }

            if (integrationPtI < 0 || integrationPtI >= statesPtr->size())
            {
                FatalErrorInFunction
                    << "Requested integration point " << integrationPtI
                    << " but valid range is [0, " << (statesPtr->size() - 1)
                    << "]." << exit(FatalError);
            }

            return (*statesPtr)[integrationPtI];
        }

        //- Return available algebraic variable names for this model
        wordList algebraicVariableNames() const
        {
            const label nAlg = ioNumAlgebraic();
            if (nAlg <= 0)
            {
                return wordList();
            }

            const char* const* names = ioAlgebraicNames();
            if (!names)
            {
                return wordList();
            }

            wordList result(nAlg);
            for (label i = 0; i < nAlg; ++i)
            {
                result[i] = word(names[i]);
            }
            return result;
        }

        //- Return algebraic values at an integration point
        scalarField algebraicVariableValues(const label integrationPtI = 0) const
        {
            const auto* algebraicPtr = ioAlgebraicPtr();
            if (!algebraicPtr || algebraicPtr->empty())
            {
                return scalarField();
            }

            if (integrationPtI < 0 || integrationPtI >= algebraicPtr->size())
            {
                FatalErrorInFunction
                    << "Requested integration point " << integrationPtI
                    << " but valid range is [0, " << (algebraicPtr->size() - 1)
                    << "]." << exit(FatalError);
            }

            return (*algebraicPtr)[integrationPtI];
        }

        // Tissue/variant metadata used by model-side setup code.
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

        virtual List<word> supportedDimensions() const
        {
            return {"1D", "2D", "3D"};
        }


//-------------------------------Print Export and write ---------------------//
        // Write and export from internal variables functions
        // singleCell
        virtual void writeHeader(OFstream& os) const
        {
            if (!hasIOMetadata())
            {
                FatalErrorInFunction
                    << "I/O metadata not available for ionic model "
                    << typeName
                    << ". Override writeHeader()/write() or provide io* hooks."
                    << exit(FatalError);
            }

            const wordList names = exportedFieldNames();

            if (!names.empty())
            {
                ionicModelIO::writeSelectedHeader(os, names);
            }
            else
            {
                ionicModelIO::writeHeader
                (
                    os,
                    ioStateNames(),
                    ioNumStates(),
                    ioAlgebraicNames(),
                    ioNumAlgebraic(),
                    fullWritePlanCache_
                );
            }
        }

        virtual void write(const scalar t, OFstream& os) const
        {
            if (!hasIOMetadata())
            {
                FatalErrorInFunction
                    << "I/O metadata not available for ionic model "
                    << typeName
                    << ". Override writeHeader()/write() or provide io* hooks."
                    << exit(FatalError);
            }

            const auto* statesPtr = ioStatesPtr();
            const auto* algebraicPtr = ioAlgebraicPtr();
            const auto* ratesPtr = ioRatesPtr();

            if (!statesPtr || !algebraicPtr || !ratesPtr)
            {
                FatalErrorInFunction
                    << "I/O storage pointers are not available for ionic model "
                    << typeName
                    << ". Override write()/exportStates() or provide io* hooks "
                    << "(including ioRatesPtr)."
                    << exit(FatalError);
            }

            const wordList names = exportedFieldNames();

            if (!names.empty())
            {
                ionicModelIO::writeSelected
                (
                    t,
                    os,
                    *statesPtr,
                    *algebraicPtr,
                    names,
                    ioStateNames(),
                    ioNumStates(),
                    ioAlgebraicNames(),
                    ioNumAlgebraic(),
                    writeSelectedPlanCache_,
                    *ratesPtr,
                    ioVmTransform()
                );
            }
            else
            {
                ionicModelIO::write
                (
                    t,
                    os,
                    *statesPtr,
                    *algebraicPtr,
                    *ratesPtr,
                    ioVmTransform()
                );
            }
        }

        // 3D vol fields(defined for states, algebraic and rates.)
        virtual void exportStates
        (
            PtrList<volScalarField>& outFields
        )
        {
            if (!hasIOMetadata())
            {
                FatalErrorInFunction
                    << "I/O metadata not available for ionic model "
                    << typeName
                    << ". Override exportStates() or provide io* hooks."
                    << exit(FatalError);
            }

            const auto* statesPtr = ioStatesPtr();
            const auto* algebraicPtr = ioAlgebraicPtr();
            const auto* ratesPtr = ioRatesPtr();

            if (!statesPtr || !algebraicPtr || !ratesPtr)
            {
                FatalErrorInFunction
                    << "I/O storage pointers are not available for ionic model "
                    << typeName
                    << ". Override exportStates() or provide io* hooks "
                    << "(including ioRatesPtr)."
                    << exit(FatalError);
            }

            ionicModelIO::exportStateFields
            (
                *statesPtr,
                *algebraicPtr,
                *ratesPtr,
                exportedFieldNames(),
                ioStateNames(),
                ioNumStates(),
                ioAlgebraicNames(),
                ioNumAlgebraic(),
                exportSelectedPlanCache_,
                outFields
            );
        }

        virtual void exportFields
        (
            const wordList& fieldNames,
            PtrList<volScalarField>& outFields
        ) const;

        virtual void importFields
        (
            const volScalarField& Vm,
            const wordList& fieldNames,
            const PtrList<volScalarField>& inFields
        );

        // Printing function
        virtual void debugPrintFields
        (
            label cellI,
            scalar t1,
            scalar t2,
            scalar step
        ) const
        {
            if (!hasIOMetadata())
            {
                FatalErrorInFunction
                    << "I/O metadata not available for ionic model "
                    << typeName
                    << ". Override debugPrintFields() or provide io* hooks."
                    << exit(FatalError);
            }

            const auto* statesPtr = ioStatesPtr();
            const auto* algebraicPtr = ioAlgebraicPtr();
            const auto* ratesPtr = ioRatesPtr();

            if (!statesPtr || !algebraicPtr || !ratesPtr)
            {
                FatalErrorInFunction
                    << "I/O storage pointers are not available for ionic model "
                    << typeName
                    << ". Override debugPrintFields() or provide io* hooks "
                    << "(including ioRatesPtr)."
                    << exit(FatalError);
            }

            const wordList names = debugPrintedNames();

            ionicModelIO::debugPrintFields
            (
                *statesPtr,
                *algebraicPtr,
                *ratesPtr,
                names,
                ioStateNames(),
                ioNumStates(),
                ioAlgebraicNames(),
                ioNumAlgebraic(),
                debugSelectedPlanCache_,
                cellI,
                t1,
                t2,
                step
            );
        }

        // Return names of variables to export
        virtual wordList exportedFieldNames() const
        {
            if (!variableExport_.size())
            {
                return wordList();
            }

            if (!hasIOMetadata())
            {
                return variableExport_;
            }

            const wordList& filtered = ionicModelIO::exportedFieldNamesRef
            (
                variableExport_,
                ioStateNames(),
                ioNumStates(),
                ioAlgebraicNames(),
                ioNumAlgebraic(),
                exportedNamesCache_
            );

            return filtered;
        }
        virtual wordList debugPrintedNames() const
        {
            if (!debugVarNames_.size())
            {
                return wordList();
            }

            if (!hasIOMetadata())
            {
                return debugVarNames_;
            }

            const wordList& filtered = ionicModelIO::exportedFieldNamesRef
            (
                debugVarNames_,
                ioStateNames(),
                ioNumStates(),
                ioAlgebraicNames(),
                ioNumAlgebraic(),
                debugNamesCache_
            );

            return filtered;
        }

        virtual void sweepCurrent
        (
            const word& currentName,
            scalar Vmin,
            scalar Vmax,
            label nPts,
            const fileName& outputFile
        ) const
        {
            FatalErrorInFunction
                << "sweepCurrent not implemented for ionic model "
                << typeName << exit(FatalError);
        }

        virtual wordList availableSweepCurrents() const
        {
            return wordList();
        }
        bool utilitiesMode() const;
        label sampleIntegrationPoint(const label nPoints) const;


    };  
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
