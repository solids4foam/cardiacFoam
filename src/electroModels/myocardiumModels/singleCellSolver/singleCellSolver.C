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

#include "singleCellSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "dimVoltage.H"
#include "stimulusIO.H"
#include "ionicModelIO.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// OverrideTypeName("singleCellSolver") in the header sets typeName_() = "singleCellSolver".
// defineTypeNameWithName registers the static member with that string; the plain
// defineTypeNameAndDebug(singleCellSolver, 0) would use #singleCellSolver and overwrite it.
defineTypeNameWithName(singleCellSolver, "singleCellSolver");
defineDebugSwitch(singleCellSolver, 0);
addToRunTimeSelectionTable(electroModel, singleCellSolver, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleCellSolver::singleCellSolver(Time& runTime, const word& region)
:
    electroModel(typeName, runTime, region),
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
    verificationModelPtr_
    (
        electroVerificationModel::New
        (
            electroProperties()
        )
    ),
    preProcessFieldNames_
    (
        verificationModelPtr_
            ? verificationModelPtr_->preProcessFieldNames(*ionicModelPtr_)
            : wordList()
    ),
    preProcessFields_(),
    postProcessFieldNames_
    (
        verificationModelPtr_
            ? verificationModelPtr_
                ->requiredPostProcessFieldNames(*ionicModelPtr_)
            : wordList()
    ),
    postProcessFields_(),
    outputPtr_(),
    activeTensionModelPtr_(),
    outputTaPtr_(),
    lambdaField_(1, 1.0),
    TaField_(1, 0.0),
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
    outFields_(),
    dummyIonicCurrentField_(1, 0.0)
{
    const fileName outputDir(runTime.path() / "postProcessing");
    mkDir(outputDir);

    const word outputSuffix =
        ionicModelIO::constantOverrideOutputSuffix(electroProperties());

    word outputName =
        ionicModelPtr_->type()
      + "_"
      + ionicModelPtr_->tissueName()
      + "_"
      + stimulusIO::protocolSuffix(electroProperties());

    if (!outputSuffix.empty())
    {
        outputName += "_";
        outputName += outputSuffix;
    }

    const fileName outFile
    (
        outputDir
      / (outputName + ".txt")
    );

    outputPtr_.reset(new OFstream(outFile));
    OFstream& output = outputPtr_.ref();

    output.setf(std::ios::fixed);
    output.precision(7);

    const wordList exportNames = ionicModelPtr_->exportedFieldNames();

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

    if (verificationModelPtr_)
    {
        electroVerificationModel::allocateFields
        (
            preProcessFieldNames_,
            mesh(),
            "preProcess_",
            preProcessFields_
        );

        electroVerificationModel::allocateFields
        (
            postProcessFieldNames_,
            mesh(),
            "postProcess_",
            postProcessFields_
        );

        verificationModelPtr_->preProcess
        (
            *ionicModelPtr_, Vm_, preProcessFields_
        );
    }
    ionicModelPtr_->writeHeader(output);

    if (electroProperties().found("activeTensionModel"))
    {
        activeTensionModelPtr_ = activeTensionModel::New
        (
            electroProperties(),
            1
        );
        activeTensionModelPtr_->setElectromechanicalSignalProvider(*ionicModelPtr_);
        activeTensionModelPtr_->validateProvider();

        const fileName outFileTa
        (
            outputDir
          / (outputName + "_Ta.txt")
        );
        outputTaPtr_.reset(new OFstream(outFileTa));
        outputTaPtr_->setf(std::ios::fixed);
        outputTaPtr_->precision(7);

        activeTensionModelPtr_->writeHeader(outputTaPtr_.ref());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool singleCellSolver::evolve()
{
    const scalar t0 = runTime().value() - runTime().deltaTValue();
    const scalar dt = runTime().deltaTValue();
    const scalar t1 = runTime().value();

    dummyIonicCurrentField_ = 0.0;
    ionicModelPtr_->solveODE
    (
        t0,
        dt,
        Vm_.internalField(),
        dummyIonicCurrentField_
    );

    if (activeTensionModelPtr_)
    {
        activeTensionModelPtr_->calculateTension(runTime().value(), dt, lambdaField_, TaField_);
    }

    const bool shouldPostProcess =
        verificationModelPtr_
      ? verificationModelPtr_->shouldPostProcess(*ionicModelPtr_, Vm_)
      : false;

    if ((runTime().outputTime() || shouldPostProcess) && !outFields_.empty())
    {
        ionicModelPtr_->exportStates(outFields_);
    }

    if (shouldPostProcess && verificationModelPtr_)
    {
        if (!postProcessFields_.empty())
        {
            ionicModelPtr_->exportFields
            (
                postProcessFieldNames_,
                postProcessFields_
            );
        }

        verificationModelPtr_->postProcess
        (
            *ionicModelPtr_,
            Vm_,
            postProcessFields_
        );
    }

    if (ionicModelIO::shouldWriteStep(t0, t1, electroProperties(), false))
    {
        ionicModelPtr_->write(runTime().value(), outputPtr_.ref());
        if (activeTensionModelPtr_)
        {
            activeTensionModelPtr_->write(runTime().value(), outputTaPtr_.ref());
        }
    }

    return true;
}


void singleCellSolver::end()
{
    runTime().printExecutionTime(Info);

    Info<< "Results written to: " << outputPtr_->name() << nl
        << "Format: [Time STATES ALGEBRAIC RATES]" << endl;

    if (activeTensionModelPtr_)
    {
        Info<< "Active tension results written to: " << outputTaPtr_->name() << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
