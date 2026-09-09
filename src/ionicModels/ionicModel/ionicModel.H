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

#include "electromechanicalSignalProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class ionicModel Declaration
\*---------------------------------------------------------------------------*/

class ionicModel
:
    public ODESystem,
    public ElectromechanicalSignalProvider
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

        //- Flag for biological sex (0=neutral, 1=male, 2=female)
        mutable label sex_;

        //- Flag for enable the update of the transmembrane voltage within the
        mutable Switch solveVmWithinODESolver_;

        //- External stimulus protocol (separate from model constants)
        StimulusProtocol stimulus_;

        //- dVm/dt (V/s) for the current step, or null to hold Vm frozen
        const scalarField* VmRatePtr_;

        //- dVm/dt of the active integration point, in model state/time units
        scalar activeVmRate_;


    // Private Member Functions

        //- No copy construct
        ionicModel(const ionicModel&) = delete;

        //- No copy assignment
        void operator=(const ionicModel&) = delete;

    protected:
        virtual void setTissueFromDict();
        void setTissue(label t) { tissue_ = t; }
        virtual void setSexFromDict();
        void setSex(label s) { sex_ = s; }
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

    public:
        // Derived models can override these and then rely on base
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

        virtual scalarField* ioMutableConstantsPtr() const
        {
            return nullptr;
        }

        virtual bool supportsRestartState() const
        {
            return false;
        }

        virtual bool readRestartState(const fvMesh&)
        {
            return false;
        }

        virtual void writeRestartState(const fvMesh&) const
        {}

        virtual void refreshRestartState(const fvMesh&)
        {}



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

        virtual void setManufacturedSourceTerm
        (
            const scalarField& sourceTerm,
            scalar chi,
            scalar Cm,
            label sourceStart = 0
        )
        {
            (void)sourceTerm;
            (void)chi;
            (void)Cm;
            (void)sourceStart;
        }

        virtual void clearManufacturedSourceTerm()
        {}

        virtual ionicModelIO::VmTransform ioVmTransform() const
        {
            return nullptr;
        }

        bool hasIOMetadata() const
        {
            return ioStateNames() && ioNumStates() > 0;
        }

    protected:

        void applyIonicConstantOverrides() const;

public:

    //- Configuration dictionary for this model instance.
    const dictionary& dict() const
    {
        return dict_;
    }

    //- Runtime type information
    TypeName("ionicModel");



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
        ionicModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );



    //- Destructor
    virtual ~ionicModel();

    // Member Functions

    bool hasSignal(CouplingSignal s) const override;
    scalar signal(label i, CouplingSignal s) const override;


        //- ODE solver
        ODESolver& odeSolver()
        {
            if (!odeSolver_)
            {
                odeSolver_.reset(ODESolver::New(*this, dict_));
            }
            return odeSolver_();
        }

        //- Step size for ODE solver
        scalarList& step()
        {
            return step_;
        }

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        ) = 0;


    // Semi-implicit (SBDF2) reaction-diffusion coupling

        //- Does solveODE honour a supplied Vm rate?
        virtual bool supportsVmRateCoupling() const
        {
            return false;
        }

        //- Can Iion be recomputed at a given Vm without advancing the state?
        virtual bool supportsIonicCurrentEvaluation() const
        {
            return false;
        }

        //- Recompute Im from the current ODE state at the supplied Vm
        virtual void evaluateIonicCurrent
        (
            const scalar t,
            const scalarField& Vm,
            scalarField& Im
        )
        {
            (void)t;
            (void)Vm;
            (void)Im;

            NotImplemented;
        }

        //- Supply dVm/dt (V/s) for the step about to be integrated
        void setVmRate(const scalarField& dVmdt)
        {
            VmRatePtr_ = &dVmdt;
        }

        //- Drop the supplied rate; Vm reverts to being held frozen
        void clearVmRate()
        {
            VmRatePtr_ = nullptr;
            activeVmRate_ = 0.0;
        }

        bool hasVmRate() const
        {
            return VmRatePtr_ != nullptr;
        }

        //- dVm/dt at integration point i in V/s, zero if none supplied
        scalar VmRateSI(const label i) const
        {
            return VmRatePtr_ ? (*VmRatePtr_)[i] : 0.0;
        }

        //- dVm/dt of the active point, in this model's state and time units
        scalar activeVmRate() const
        {
            return activeVmRate_;
        }

    protected:

        //- Latch the Vm rate for integration point i under the given scalings
        void setActiveVmRate
        (
            const label i,
            const scalar stateScale = 1.0,
            const scalar timeScale = 1.0
        )
        {
            activeVmRate_ = VmRateSI(i)*stateScale/timeScale;
        }

    public:

        virtual void configureIonicHeterogeneity
        (
            const scalarField& transmuralDistance,
            const dictionary& heterogeneityDict
        );

        virtual void configureGradientAxisHeterogeneity
        (
            const scalarField& fieldValues,
            const dictionary& dict
        );

        //- Configure region-based heterogeneity through
        //  ionicHeterogeneityOrchestrator.
        void configureTransmuralBandHeterogeneity
        (
            const scalarField& transmuralDistance,
            const dictionary& heterogeneityDict,
            PtrList<scalarField>& heterogeneousConstants,
            PtrList<scalarField>* heterogeneousInitialStates = nullptr
        ) const;

        //- Configure one gradientAxes axis's heterogeneity.
        void configureGradientAxisHeterogeneityImpl
        (
            const scalarField& fieldValues,
            const dictionary& dict,
            PtrList<scalarField>& heterogeneousConstants
        ) const;

        virtual scalarField constantsForTissue(const label tissueFlag) const;

        virtual scalarField initialStatesForTissue(const label tissueFlag) const;

        //- Number of ODE's to solve
        virtual label nEqns() const override = 0;

        //- Derivatives of each ODE with respect to the primary unknown in that
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



        //- Getter for tissue flag
        label tissue() const
        {
            return tissue_;
        }

        //- Getter for sex flag
        label sex() const
        {
            return sex_;
        }

        virtual word sexName() const
        {
            const label s = sex_;
            return (s == 1) ? "male" :
                   (s == 2) ? "female" :
                              "neutral";
        }

        //- Return the resting membrane potential in SI units (V), derived
        virtual scalar vmRest() const
        {
            const auto* statesPtr = ioStatesPtr();
            if (!statesPtr || statesPtr->empty() || (*statesPtr)[0].empty())
            {
                return -0.084;
            }
            const scalarField& s0 = (*statesPtr)[0];
            const ionicModelIO::VmTransform transform = ioVmTransform();
            const scalar vmMv = transform ? transform(s0) : s0[0];
            return vmMv * 1e-3;   // mV → V
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
                if (!names[i] || names[i][0] == '\0')
                {
                    FatalErrorInFunction
                        << "Invalid constant metadata for ionic model "
                        << type() << ": constant index " << i
                        << " has no name. Fix the model's CONSTANTS_NAMES "
                        << "array so every NUM_CONSTANTS entry has a real "
                        << "semantic name."
                        << exit(FatalError);
                }

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
        scalarField algebraicVariableValues
        (const label integrationPtI = 0) const
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

        virtual List<word> supportedTissueTypes() const
        {
            return {"epicardialCells", "mCells", "endocardialCells", "myocyte"};
        }

        virtual List<word> supportedSexTypes() const
        {
            return {};   // empty = model does not support sex selection
        }

        //- Optional verification family tag for manufactured/reference tests.
        virtual word verificationFamily() const
        {
            return word();
        }

        //- Optional geometric dimension used by manufactured/reference tests.
        virtual label geometricDimension() const
        {
            return -1;
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
                    ioStateNames(),
                    ioNumStates(),
                    ioAlgebraicNames(),
                    ioNumAlgebraic(),
                    fullWritePlanCache_,
                    *ratesPtr,
                    ioVmTransform()
                );
            }
        }

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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
