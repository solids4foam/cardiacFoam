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

#include "Grandi.H"
#include "Grandi_2011.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelFamilyInfo.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"

#include <math.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Grandi, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, Grandi, dictionary
    );

    const ionicModelFamilyInfo& GrandiFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            GrandiCONSTANTS_NAMES,
            GrandiSTATES_NAMES,
            GrandiALGEBRAIC_NAMES,
            membrane_V,
            1000.0,
            1000.0,
            0.0,
            nullptr
        };
        return info;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Grandi::Grandi
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{

    // 🔑 First, set tissue using base logic + overrides
    ionicModel::setTissueFromDict();
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        GrandiinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue(),dict
        );
        if (!utilitiesMode())
        {
            setStimulusProtocolFromDict(dict);
        }
    }

    applyIonicConstantOverrides();

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Grandi::~Grandi()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::Grandi::supportedTissueTypes() const
{
    return {"myocyte"};
}


// ------------------------------------------------------------------------- //
//  Solve ODE with mixed singleCell implementation and 1D-3D condition
// ------------------------------------------------------------------------- //
void Foam::Grandi::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const scalar tStart = stepStartTime * 1000;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000;
    const label sampleCell = sampleIntegrationPoint(STATES_.size());

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& step = ionicModel::step()[integrationPtI];

        if (!solveVmWithinODESolver())
        {
            STATESI[membrane_V] = Vm[integrationPtI]*1000.0;
        }

        step = min(step, deltaT * 1000.0);
        odeSolver().solve(tStart, tEnd, STATESI, step);

        ::GrandicomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );
        if (integrationPtI == sampleCell)
        {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;

    }
}

void Foam::Grandi::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::GrandicomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
        tissue(),
        solveVmWithinODESolver()
    ,
            stimulusProtocol()
        );
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

const char* const* Foam::Grandi::ioStateNames() const
{
    return GrandiSTATES_NAMES;
}

const char* const* Foam::Grandi::ioConstantNames() const
{
    return GrandiCONSTANTS_NAMES;
}

const char* const* Foam::Grandi::ioAlgebraicNames() const
{
    return GrandiALGEBRAIC_NAMES;
}

void Foam::Grandi::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = GrandiDependencyMap();

    if (!depMap.found(currentName))
    {
        FatalErrorInFunction
            << "Unknown current: " << currentName << nl
            << "Available currents: " << depMap.toc() << nl
            << exit(FatalError);
    }

    const wordList& deps = depMap[currentName];
    OFstream os(outputFile);
    ionicModelIO::writeSweepHeader(os, deps);

    scalarField STATESI = STATES_[0];
    scalarField RATESI(NUM_STATES, 0.0);
    scalarField ALGI(NUM_ALGEBRAIC, 0.0);
    ionicModelIO::SelectedMapCache sweepPlanCache;

    for (label i = 0; i < nPts; ++i)
    {
        scalar V = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);

        STATESI = STATES_[0];

        STATESI[0] = V;

        ::GrandicomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps,STATESI,ALGI,
            GrandiSTATES_NAMES, NUM_STATES,
            GrandiALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::Grandi::availableSweepCurrents() const
{
    return GrandiDependencyMap().toc();
}
