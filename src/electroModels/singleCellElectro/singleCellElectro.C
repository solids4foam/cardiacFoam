/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "singleCellElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "stimulusIO.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(singleCellElectro, 0);
addToRunTimeSelectionTable(electroModel, singleCellElectro, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleCellElectro::singleCellElectro
(
    Time& runTime,
    const word& region
)
:
    electroModel(typeName, runTime, region),
    // cardiacProperties_(electroProperties().subDict("cardiacProperties")),
    cardiacProperties_(electroProperties()),
    ionicModelPtr_
    (
        ionicModel::New
        (
            electroProperties(),
            1, // one integration point
            runTime.deltaTValue(),
            true // solve Vm equation within ODE system
        )
    ),
    outputPtr_()
{
    // Create the output file
    const fileName outFile =
        runTime.path()/ionicModelPtr_->type() + "_"
      + ionicModelPtr_->tissueName() + "_"
      + stimulusIO::protocolSuffix(electroProperties()) +".txt";

    outputPtr_.reset(new OFstream(outFile));
    OFstream& output = outputPtr_.ref();

    output.setf(std::ios::fixed);
    output.precision(7);

    // Extract the names of the fields to be exported
    const wordList exportNames = ionicModelPtr_->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }

     ionicModelPtr_->writeHeader(output);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool singleCellElectro::evolve()
{
    // Create dummy Vm, Iion and state fields
    Field<Field<scalar>> dummyStates(1);
    dummyStates[0].setSize(ionicModelPtr_->nEqns());
    const scalarField dummyVmField(1, 0);
    scalarField dummyIonicCurrentField(1, 0.0);

    // Old time
    const scalar t0 = runTime().value() - runTime().deltaTValue();

    // Time step
    const scalar dt = runTime().deltaTValue();

    // Current time
    const scalar t1 = runTime().value();

    // Solve the ionic model from t0 to t1
    ionicModelPtr_->solveODE
    (
        t0,
        dt,
        dummyVmField,
        dummyIonicCurrentField,
        dummyStates
    );

    if (stimulusIO::shouldWriteStep(t0, t1, electroProperties(), false))
    {
        ionicModelPtr_->write(runTime().value(), outputPtr_.ref());
    }

    return true;
}


void singleCellElectro::end()
{
    runTime().printExecutionTime(Info);

    Info<< "Results written to: " << outputPtr_->name() << nl
        << "Format: [Time STATES ALGEBRAIC RATES]" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
