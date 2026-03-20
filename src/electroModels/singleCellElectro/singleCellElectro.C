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

#include "singleCellElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "stimulusIO.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(singleCellElectro, 0);
addToRunTimeSelectionTable(electroModel, singleCellElectro, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleCellElectro::singleCellElectro(Time& runTime, const word& region)
:
    electroModel(typeName, runTime, region),
    cardiacProperties_(electroProperties()),
    ionicModelPtr_
    (
        ionicModel::New
        (
            ionicProperties(),
            1,
            runTime.deltaTValue(),
            true
        )
    ),
    prePostProcessor_
    (
        electroModelsPrePostProcessor::New
        (
            ionicModelPtr_->type(),
            ionicProperties()
        )
    ),
    preProcessFieldNames_(prePostProcessor_->preProcessFieldNames(*ionicModelPtr_)),
    preProcessFields_(),
    outputPtr_(),
    Vm_
    (
        IOobject
        (
            "Vm",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("Vm", dimVoltage, -80.0),
        "zeroGradient"
    ),
    outFields_()
{
    const fileName outputDir(runTime.path() / "postProcessing");
    mkDir(outputDir);

    const fileName outFile
    (
        outputDir
      / (
            ionicModelPtr_->type()
          + "_"
          + ionicModelPtr_->tissueName()
          + "_"
          + stimulusIO::protocolSuffix(electroProperties())
          + ".txt"
        )
    );

    outputPtr_.reset(new OFstream(outFile));
    OFstream& output = outputPtr_.ref();

    output.setf(std::ios::fixed);
    output.precision(7);

    const wordList exportNames = ionicModelPtr_->exportedFieldNames();
    const wordList requiredPostProcessNames =
        prePostProcessor_->requiredPostProcessFieldNames(*ionicModelPtr_);

    if (!exportNames.empty())
    {
        Info<< "Exporting fields: " << exportNames << nl;
    }

    outFields_.setSize(exportNames.size());
    forAll(exportNames, i)
    {
        outFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    exportNames[i],
                    runTime.timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimless,
                "zeroGradient"
            )
        );
    }

    electroModelsPrePostProcessor::validatePostProcessFields
    (
        ionicModelPtr_->type(),
        exportNames,
        requiredPostProcessNames
    );

    electroModelsPrePostProcessor::allocateFields
    (
        preProcessFieldNames_,
        mesh(),
        "preProcess_",
        preProcessFields_
    );

    prePostProcessor_->preProcess(*ionicModelPtr_, Vm_, preProcessFields_);
    ionicModelPtr_->writeHeader(output);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool singleCellElectro::evolve()
{
    scalarField dummyIonicCurrentField(1, 0.0);

    const scalar t0 = runTime().value() - runTime().deltaTValue();
    const scalar dt = runTime().deltaTValue();
    const scalar t1 = runTime().value();

    ionicModelPtr_->solveODE(t0, dt, Vm_.internalField(), dummyIonicCurrentField);

    const bool shouldPostProcess =
        prePostProcessor_->shouldPostProcess(*ionicModelPtr_, Vm_);

    if ((runTime().outputTime() || shouldPostProcess) && !outFields_.empty())
    {
        ionicModelPtr_->exportStates(outFields_);
    }

    if (shouldPostProcess)
    {
        prePostProcessor_->postProcess(*ionicModelPtr_, Vm_, outFields_);
    }

    if (ionicModelIO::shouldWriteStep(t0, t1, electroProperties(), false))
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
