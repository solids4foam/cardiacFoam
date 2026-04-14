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

\*---------------------------------------------------------------------------*/

#ifndef activeTensionModel_H
#define activeTensionModel_H

#include "activeTensionIO.H"
#include "electromechanicalSignalProvider.H"

#include "OFstream.H"
#include "PtrList.H"
#include "Switch.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "error.H"
#include "label.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFields.H"
#include "wordList.H"

namespace Foam
{

class activeTensionModel
{
protected:

    const dictionary& dict_;
    const label nIntegrationPoints_;

    // Optional pointer to coupling signal provider (typically ionicModel or adapter)
    const ElectromechanicalSignalProvider* providerPtr_;

    //- Drive signal value for the current integration point.
    //  Set by the base calculateTension() loop before each solveAtPoint() call
    //  so th at derivatives() can read it during ODE sub-steps.
    scalar currentDriveSignal_;

    //- Time at start of current timestep (set by base before loop)
    scalar currentT_;

    //- Timestep size (set by base before loop)
    scalar currentDt_;

    wordList debugVarNames_;
    wordList exportVarNames_;

    // I/O plan caches (mutable: populated lazily on first use)
    mutable activeTensionIO::FullPlanCache    fullWritePlanCache_;
    mutable activeTensionIO::SelectedMapCache writeSelectedPlanCache_;
    mutable activeTensionIO::SelectedMapCache exportSelectedPlanCache_;
    mutable activeTensionIO::SelectedMapCache debugSelectedPlanCache_;
    mutable activeTensionIO::ExportedNamesCache exportedNamesCache_;
    mutable activeTensionIO::ExportedNamesCache debugNamesCache_;

    // ---- Optional metadata hooks for generic I/O ----
    // Derived classes override these to enable base-class write/export/debug.

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

    bool hasIOMetadata() const
    {
        return ioStateNames() != nullptr && ioNumStates() > 0;
    }


public:

    TypeName("activeTensionModel");

    // Runtime selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        activeTensionModel,
        dictionary,
        (const dictionary& dict, const label nIntegrationPoints),
        (dict, nIntegrationPoints)
    );

    // Construct from dictionary
    activeTensionModel(const dictionary& dict, const label nIntegrationPoints)
    :
        dict_(dict),
        nIntegrationPoints_(nIntegrationPoints),
        providerPtr_(nullptr),
        currentDriveSignal_(0.0),
        currentT_(0.0),
        currentDt_(0.0)
    {
        // Required schema:
        // outputVariables
        // {
        //   activeTension
        //   {
        //     export (...);
        //     debug  (...);
        //   }
        // }
        if (dict_.found("outputVariables"))
        {
            const dictionary& outDict = dict_.subDict("outputVariables");
            if (outDict.found("activeTension"))
            {
                const dictionary& atOut = outDict.subDict("activeTension");

                if (atOut.found("export"))
                {
                    atOut.lookup("export") >> exportVarNames_;
                }

                if (atOut.found("debug"))
                {
                    atOut.lookup("debug") >> debugVarNames_;
                }
            }
        }
    }

    virtual ~activeTensionModel() = default;

    // Factory selector
    static autoPtr<activeTensionModel> New
    (
        const dictionary& dict,
        const label nIntegrationPoints
    );

    // Requirements declared by derived models
    struct Requirements
    {
        bool needVm  = false;
        bool needCai = false;
    };

    // Derived models declare their required inputs
    virtual Requirements requirements() const = 0;

    // Attach provider (call once from solver/harness)
    virtual void setElectromechanicalSignalProvider
    (
        const ElectromechanicalSignalProvider& p
    )
    {
        providerPtr_ = &p;
    }

    // Validate that the provider satisfies requirements (call once after wiring)
    virtual void validateProvider() const;

    // Convenience accessors for derived classes (fail-fast + clear errors)
    inline const ElectromechanicalSignalProvider& provider() const
    {
        if (!providerPtr_)
        {
            FatalErrorInFunction
                << "ElectromechanicalSignalProvider not set. "
                << "Call setElectromechanicalSignalProvider(...) before correct()."
                << abort(FatalError);
        }
        return *providerPtr_;
    }

    label nIntegrationPoints() const
    {
        return nIntegrationPoints_;
    }

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

    scalarField constantVariableValues() const
    {
        const auto* constantsPtr = ioConstantsPtr();
        if (!constantsPtr)
        {
            return scalarField();
        }

        return *constantsPtr;
    }

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

    scalarField algebraicVariableValues(const label integrationPtI = 0) const
    {
        const auto* algebraicPtr = ioAlgebraicPtr();
        if (!algebraicPtr || algebraicPtr->empty())
        {
            return scalarField();
        }

        if
        (
            integrationPtI < 0
         || integrationPtI >= algebraicPtr->size()
        )
        {
            FatalErrorInFunction
                << "Requested integration point " << integrationPtI
                << " but valid range is [0, "
                << (algebraicPtr->size() - 1) << "]."
                << exit(FatalError);
        }

        return (*algebraicPtr)[integrationPtI];
    }

    wordList rateVariableNames() const
    {
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
            result[i] = word("RATES_") + word(names[i]);
        }
        return result;
    }

    scalarField rateVariableValues(const label integrationPtI = 0) const
    {
        const auto* ratesPtr = ioRatesPtr();
        if (!ratesPtr || ratesPtr->empty())
        {
            return scalarField();
        }

        if (integrationPtI < 0 || integrationPtI >= ratesPtr->size())
        {
            FatalErrorInFunction
                << "Requested integration point " << integrationPtI
                << " but valid range is [0, " << (ratesPtr->size() - 1)
                << "]." << exit(FatalError);
        }

        return (*ratesPtr)[integrationPtI];
    }

    // Return names of variables to export (filtered against available names)
    virtual wordList exportedFieldNames() const
    {
        if (!exportVarNames_.size())
            return wordList();

        if (!hasIOMetadata())
            return exportVarNames_;

        return activeTensionIO::exportedFieldNamesRef
        (
            exportVarNames_,
            ioStateNames(), ioNumStates(),
            ioAlgebraicNames(), ioNumAlgebraic(),
            exportedNamesCache_
        );
    }

    virtual wordList debugPrintedNames() const
    {
        if (!debugVarNames_.size())
            return wordList();

        if (!hasIOMetadata())
            return debugVarNames_;

        return activeTensionIO::exportedFieldNamesRef
        (
            debugVarNames_,
            ioStateNames(), ioNumStates(),
            ioAlgebraicNames(), ioNumAlgebraic(),
            debugNamesCache_
        );
    }

    // Flat list of all available variable names (states + algebraic + rates).
    // Implemented in base using io* hooks; derived classes may override.
    virtual wordList availableFieldNames() const;

    // Flat field values at the given integration point.
    // Implemented in base using io* hooks; derived classes may override.
    virtual void fieldValues(const label integrationPtI, scalarField& values) const;

    //- Return the CouplingSignal this model reads as its drive input.
    //  Active-tension feedback accepts Vm or Cai.
    virtual CouplingSignal driveSignal() const = 0;

    //- Advance ODE and write active tension for a single integration point.
    //  Called by the base calculateTension() loop.
    //  Derived classes must NOT call provider() directly — driveVal is already resolved.
    //
    //  i        : integration-point (cell) index
    //  driveVal : coupling signal value at this cell (Vm / Cai)
    //  lambda   : stretch at this cell (ignored by models that don't use it)
    //  Ta       : output — active tension at this cell
    virtual void solveAtPoint
    (
        label i,
        scalar driveVal,
        scalar lambda,
        scalar& Ta
    ) = 0;

    //- Base implementation: queries provider, loops over cells, calls solveAtPoint().
    //  Derived classes must NOT override this.
    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    );

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const = 0;

    // Write header for single-cell CSV output
    virtual void writeHeader(OFstream& os) const;

    // Write one data row for single-cell CSV output (integration point 0)
    virtual void write(const scalar t, OFstream& os) const;

    // Export state fields to volScalarField output list
    virtual void exportStates(PtrList<volScalarField>& outFields);

    // Print debug fields for a given integration point
    virtual void debugPrintFields
    (
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    ) const;
};

} // End namespace Foam

#endif
