/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Base class for active tension models
\*---------------------------------------------------------------------------*/

#ifndef activeTensionModel_H
#define activeTensionModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "PtrList.H"
#include "scalarField.H"
#include "label.H"
#include "Switch.H"
#include "error.H"
#include "wordList.H"
#include "OFstream.H"
#include "volFields.H"

#include "activeTensionIO.H"
#include "electromechanicsFeedbackProvider.H"

namespace Foam
{

class activeTensionModel
{
protected:

    const dictionary& dict_;
    const label nIntegrationPoints_;

    // Optional pointer to coupling signal provider (typically ionicModel or adapter)
    const CouplingSignalProvider* providerPtr_;

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
        providerPtr_(nullptr)
    {
        if (dict_.found("ActiveTensionModel_debugPrint"))
            dict_.lookup("ActiveTensionModel_debugPrint") >> debugVarNames_;
        else if (dict_.found("debugPrintVariables"))
            dict_.lookup("debugPrintVariables") >> debugVarNames_;

        if (dict_.found("ActiveTensionModel_export"))
            dict_.lookup("ActiveTensionModel_export") >> exportVarNames_;
        else if (dict_.found("exportedVariables"))
            dict_.lookup("exportedVariables") >> exportVarNames_;
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
        bool needAct = false;
        bool needCai = false;
    };

    // Derived models declare their required inputs
    virtual Requirements requirements() const = 0;

    // Attach provider (call once from solver/harness)
    virtual void setCouplingSignalProvider(const CouplingSignalProvider& p)
    {
        providerPtr_ = &p;
    }

    // Validate that the provider satisfies requirements (call once after wiring)
    virtual void validateProvider() const;

    // Convenience accessors for derived classes (fail-fast + clear errors)
    inline const CouplingSignalProvider& provider() const
    {
        if (!providerPtr_)
        {
            FatalErrorInFunction
                << "CouplingSignalProvider not set. "
                << "Call setCouplingSignalProvider(...) before correct()."
                << abort(FatalError);
        }
        return *providerPtr_;
    }

    label nIntegrationPoints() const
    {
        return nIntegrationPoints_;
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

    // Main API: update active tension for all points
    virtual void calculateTension
    (
        const scalar t,
        const scalar dt,
        const scalarField& lambda,
        scalarField& Ta
    ) = 0;

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
