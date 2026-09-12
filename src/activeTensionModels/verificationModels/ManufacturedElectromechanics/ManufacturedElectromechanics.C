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

#include "ManufacturedElectromechanics.H"
#include "ManufacturedElectromechanics_2026.H"

#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "word.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

defineTypeNameAndDebug(ManufacturedElectromechanics, 0);

addToRunTimeSelectionTable
(
    activeTensionModel,
    ManufacturedElectromechanics,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const char* const* ManufacturedElectromechanics::ioConstantNames() const
{
    return ManufacturedElectromechanicsCONSTANTS_NAMES;
}


const char* const* ManufacturedElectromechanics::ioStateNames() const
{
    return ManufacturedElectromechanicsSTATES_NAMES;
}


const char* const* ManufacturedElectromechanics::ioAlgebraicNames() const
{
    return ManufacturedElectromechanicsALGEBRAIC_NAMES;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

ManufacturedElectromechanics::ManufacturedElectromechanics
(
    const dictionary& dict,
    const label nIntegrationPoints
)
:
    activeTensionModel(dict, nIntegrationPoints),
    STATES_(nIntegrationPoints),
    ALGEBRAIC_(nIntegrationPoints),
    RATES_(nIntegrationPoints),
    CONSTANTS_(NUM_CONSTANTS, 0.0)
{
    scalarField protoStates(NUM_STATES, 0.0);
    scalarField protoRates(NUM_STATES, 0.0);

    ManufacturedElectromechanics2026initConsts
    (
        CONSTANTS_.data(),
        protoRates.data(),
        protoStates.data()
    );

    if (dict.found("constants"))
    {
        const dictionary& cDict = dict.subDict("constants");
        for (label k = 0; k < NUM_CONSTANTS; ++k)
        {
            const word name(ManufacturedElectromechanicsCONSTANTS_NAMES[k]);
            if (cDict.found(name))
            {
                CONSTANTS_[k] = cDict.get<scalar>(name);
            }
        }
    }

    if (CONSTANTS_[MMS_AC_V0] <= SMALL)
    {
        FatalErrorInFunction
            << "ManufacturedElectromechanics constant V0 must be positive. "
            << "Current value: " << CONSTANTS_[MMS_AC_V0]
            << abort(FatalError);
    }

    if (dict.found("initialStates"))
    {
        const dictionary& sDict = dict.subDict("initialStates");
        for (label k = 0; k < NUM_STATES; ++k)
        {
            const word name(ManufacturedElectromechanicsSTATES_NAMES[k]);
            if (sDict.found(name))
            {
                protoStates[k] = sDict.get<scalar>(name);
            }
        }
    }

    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(protoStates));
        ALGEBRAIC_.set(i, new scalarField(NUM_ALGEBRAIC, 0.0));
        RATES_.set(i, new scalarField(protoRates));
    }

    Info<< nl
        << "Initialize ManufacturedElectromechanics constants:" << nl
        << CONSTANTS_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

activeTensionModel::Requirements
ManufacturedElectromechanics::requirements() const
{
    Requirements req;
    req.needVm = true;
    return req;
}


void ManufacturedElectromechanics::solveAtPoint
(
    const label i,
    const scalar driveVal,
    const scalar lambda,
    scalar& Ta
)
{
    scalarField& STATESI = STATES_[i];
    scalarField& ALGEBRAICI = ALGEBRAIC_[i];
    scalarField& RATESI = RATES_[i];

    ALGEBRAICI[MMS_AV_Vm] = driveVal;
    ALGEBRAICI[MMS_AV_lambda] = lambda;

    ManufacturedElectromechanics2026computeVariables
    (
        currentT_ + currentDt_,
        CONSTANTS_.data(),
        RATESI.data(),
        STATESI.data(),
        ALGEBRAICI.data()
    );

    Ta = STATESI[MMS_STATE_Ta];

    if (!debugPrintedNames().empty() && i == 0)
    {
        debugPrintFields(i, currentT_, currentT_ + currentDt_, currentDt_);
    }
}


void ManufacturedElectromechanics::derivatives
(
    const scalar,
    const scalarField&,
    scalarField& dydt
) const
{
    dydt.setSize(NUM_STATES);
    dydt[MMS_STATE_Ta] = 0.0;
}

} // End namespace Foam

// ************************************************************************* //
