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
#include "scalarField.H"
#include "label.H"
#include "Switch.H"
#include "error.H"
#include "wordList.H"

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
    bool needVm       = false;
    bool needAct      = false;
    bool needCai      = false;
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

    virtual wordList debugPrintedNames() const
    {
        return debugVarNames_.size()
            ? debugVarNames_
            : wordList();
    }

    virtual wordList exportedFieldNames() const
    {
        return exportVarNames_.size()
            ? exportVarNames_
            : wordList();
    }

    // Main API: update active tension for all points
    // lambda: stretch for each point (use all-ones for EP-only tests)
    // Ta: output active tension for each point (size must be nIntegrationPoints)
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

     

};

} // End namespace Foam

#endif
