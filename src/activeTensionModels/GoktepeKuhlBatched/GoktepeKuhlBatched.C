/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

\*---------------------------------------------------------------------------*/

#include "GoktepeKuhlBatched.H"
#include "GoktepeKuhlBatch.H"
#include "addToRunTimeSelectionTable.H"
#include "../GoktepeKuhl/GoktepeKuhl_2004Names.H"

namespace Foam
{
    defineTypeNameAndDebug(GoktepeKuhlBatched, 0);
    addToRunTimeSelectionTable(activeTensionModel, GoktepeKuhlBatched, dictionary);
}


// * * * * * * * * * * * * * * * * io* hooks  * * * * * * * * * * * * * * * //

const char* const* Foam::GoktepeKuhlBatched::ioStateNames() const
{
    // Need to cast to match const char* const*
    static const char* const names[NUM_STATES] = { "Ta" };
    return names;
}

const char* const* Foam::GoktepeKuhlBatched::ioConstantNames() const
{
    static const char* const names[NUM_CONSTANTS] = {
        "AC_Vr", "AC_eInfty", "AC_e0", "AC_eXi", "AC_Vshift", "AC_kTa"
    };
    return names;
}

const char* const* Foam::GoktepeKuhlBatched::ioAlgebraicNames() const
{
    static const char* const names[NUM_ALGEBRAIC] = {
        "AV_e", "AV_Vm", "AV_u"
    };
    return names;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GoktepeKuhlBatched::GoktepeKuhlBatched
(
    const dictionary& dict,
    const label num
)
:
    batchedActiveTensionModel(dict, num, NUM_STATES, NUM_ALGEBRAIC),
    CONSTANTS_(NUM_CONSTANTS, 0.0)
{
    const word requestedSignal = dict_.lookupOrDefault<word>("couplingSignal", "Vm");

    if (!(requestedSignal == "Vm" || requestedSignal == "vm"))
    {
        FatalErrorInFunction
            << "Unknown GoktepeKuhlBatched 'couplingSignal' value: "
            << requestedSignal << nl
            << "Valid option is: Vm."
            << abort(FatalError);
    }

    Info<< nl << "Initialize GoktepeKuhlBatched constants:" << nl;

    // Default Initialization
    CONSTANTS_[AC_Vr] = -80.0;
    CONSTANTS_[AC_eInfty] = 10.0;
    CONSTANTS_[AC_e0] = 1.0;
    CONSTANTS_[AC_eXi] = 1.0;
    CONSTANTS_[AC_Vshift] = 0.0;
    CONSTANTS_[AC_kTa] = 47.9;

    scalarField protoStates(NUM_STATES, 0.0);

    if (dict.found("constants"))
    {
        const dictionary& cDict = dict.subDict("constants");
        for (label k = 0; k < NUM_CONSTANTS; ++k)
        {
            const word name(ioConstantNames()[k]);
            if (cDict.found(name)) CONSTANTS_[k] = cDict.get<scalar>(name);
        }
    }

    if (dict.found("initialStates"))
    {
        const dictionary& sDict = dict.subDict("initialStates");
        for (label k = 0; k < NUM_STATES; ++k)
        {
            const word name(ioStateNames()[k]);
            if (sDict.found(name)) protoStates[k] = sDict.get<scalar>(name);
        }
    }

    for (label cellI = 0; cellI < nCells_; ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = protoStates[stateI];
        }
    }

    Info<< CONSTANTS_ << nl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GoktepeKuhlBatched::evaluateHotPathStateForCell
(
    const label cellI,
    const scalar modelTime,
    const scalar driveSignal,
    const scalar lambda,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    GoktepeKuhlComputeVariablesBatch
    (
        driveSignal,
        CONSTANTS_.data(),
        1,
        0,
        const_cast<scalar*>(stateValues.begin()),
        rateValues.begin(),
        algebraicValues.begin()
    );
}


bool Foam::GoktepeKuhlBatched::rushLarsenParametersForCell
(
    const label cellI,
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& algebraicValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI == ::Ta)
    {
        const scalar eFactor = algebraicValues[::AV_e];
        if (eFactor > VSMALL)
        {
            tau = 1.0 / eFactor;
            steadyState = CONSTANTS_[AC_kTa] * algebraicValues[::AV_u];
            return true;
        }
    }
    return false;
}

// ************************************************************************* //
