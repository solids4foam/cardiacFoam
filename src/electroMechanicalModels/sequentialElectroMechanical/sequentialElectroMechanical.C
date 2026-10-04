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

#include "sequentialElectroMechanical.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroMechanicalModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sequentialElectroMechanical, 0);
addToRunTimeSelectionTable
(
    electroMechanicalModel, sequentialElectroMechanical, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sequentialElectroMechanical::sequentialElectroMechanical
(
    Time& runTime,
    const word& region
)
:
    electroMechanicalModel(typeName, runTime, region),
    Ta_
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solid().mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        "zeroGradient"
    ),
    TaScale_
    (
        electroMechanicalProperties().lookupOrDefault<scalar>("TaScale", 1e3)
    ),
    lambdaField_(electro().mesh().nCells(), 1.0),
    activeTensionModel_
    (
        activeTensionModel::New
        (
            electroMechanicalProperties(),
            electro().mesh().nCells()
        )
    ),
    verificationModelPtr_(),
    activeTensionRequirements_(activeTensionModel_->requirements())
{
    const ElectromechanicalSignalProvider* prov = electro().provider();

    if (prov)
    {
        activeTensionModel_->setElectromechanicalSignalProvider(*prov);
    }

    activeTensionModel_->validateProvider();

    if (solid().mesh().nCells() != electro().mesh().nCells())
    {
        FatalErrorInFunction
            << "sequentialElectroMechanical requires conforming meshes "
            << "(same cell count). Solid has "
            << solid().mesh().nCells() << " cells, electro has "
            << electro().mesh().nCells() << " cells."
            << abort(FatalError);
    }

    if (!solid().mesh().foundObject<volVectorField>("f0"))
    {
        new volVectorField
        (
            IOobject
            (
                "f0",
                runTime.timeName(),
                solid().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            solid().mesh()
        );

        Info<< "    Registered f0 in solid objectRegistry." << nl << endl;
    }

    if (activeTensionRequirements_.needsLambda)
    {
        if (!solid().mesh().foundObject<volVectorField>("D"))
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field D "
                << "is not in the solid objectRegistry."
                << abort(FatalError);
        }
        if (!solid().mesh().foundObject<volVectorField>("f0"))
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field f0 "
                << "is not in the solid objectRegistry."
                << abort(FatalError);
        }
    }

    const bool activeTensionRestarted =
        activeTensionModel_->readRestartState(solid().mesh());

    if (activeTensionRestarted)
    {
        activeTensionModel_->refreshRestartState(solid().mesh());
        scalarField restartTa(Ta_.primitiveField());
        if (activeTensionModel_->restartTension(restartTa))
        {
            if (TaScale_ != 1.0)
            {
                restartTa *= TaScale_;
            }
            Ta_.primitiveFieldRef() = restartTa;
            Ta_.correctBoundaryConditions();
        }
    }

    const myocardiumPrePacing* prePacingPtr =
        electro().mesh().foundObject<myocardiumPrePacing>
        (
            myocardiumPrePacing::typeName
        )
      ? &electro().mesh().lookupObject<myocardiumPrePacing>
        (
            myocardiumPrePacing::typeName
        )
      : nullptr;

    bool allRegionsPrePaced = false;
    if (prePacingPtr)
    {
        allRegionsPrePaced = true;
        forAll(prePacingPtr->regionNames(), regionI)
        {
            allRegionsPrePaced =
                allRegionsPrePaced && prePacingPtr->paced(regionI);
        }
    }

    if
    (
        prov
     && activeTensionRequirements_.needCai
     && !activeTensionRestarted
     && !allRegionsPrePaced
    )
    {
        scalarField restingCai(electro().mesh().nCells());
        forAll(restingCai, cellI)
        {
            restingCai[cellI] = prov->signal(cellI, CouplingSignal::CAI);
        }
        activeTensionModel_->preconditionToRestingState(restingCai);
    }

    if (prov && prePacingPtr && !activeTensionRestarted)
    {
        prePaceActiveTension(*prePacingPtr);
    }

    if
    (
        electromechanicalVerificationModel::configured
        (
            electroMechanicalProperties()
        )
    )
    {
        verificationModelPtr_ =
            electromechanicalVerificationModel::New
            (
                electroMechanicalProperties()
            );
    }

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->initialize
        (
            const_cast<volScalarField&>(electro().Vm()),
            solid().D()
        );
    }

    Info<< "    Active tension model: "
        << activeTensionModel_->type() << nl
        << "    TaScale (model units -> Pa): " << TaScale_ << nl
        << "    Integration points: " << electro().mesh().nCells() << nl
        << endl;
}


void sequentialElectroMechanical::writeFields(const Time& runTime)
{
    electroMechanicalModel::writeFields(runTime);
    electro().writeRestartState();
    activeTensionModel_->writeRestartState(solid().mesh());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sequentialElectroMechanical::prePaceActiveTension
(
    const myocardiumPrePacing& prePacing
)
{
    // Single tension cell at lambda = 1, driven by its ionic cell's signal
    class isometricTension
    :
        public prePacingCoupling
    {
        activeTensionModel& model_;
        const scalarField lambda_;
        scalarField Ta_;

    public:

        explicit isometricTension(activeTensionModel& model)
        :
            model_(model),
            lambda_(1, 1.0),
            Ta_(1, 0.0)
        {}

        virtual void advance(const scalar t, const scalar dt)
        {
            model_.calculateTension(t, dt, lambda_, Ta_);
        }

        virtual scalarField state() const
        {
            return (*model_.states())[0];
        }

        scalar Ta() const
        {
            return Ta_[0];
        }
    };

    const labelList& cellRegion = prePacing.cellRegion();
    const PtrList<scalarField>* tissueStatesPtr = activeTensionModel_->states();

    if
    (
        !tissueStatesPtr
     || tissueStatesPtr->size() != cellRegion.size()
    )
    {
        FatalErrorInFunction
            << "Active tension model '" << activeTensionModel_->type()
            << "' does not expose states for all " << cellRegion.size()
            << " myocardium cells; prePacing cannot seed it. Remove "
            << "constant/prePacingProperties."
            << exit(FatalError);
    }

    const wordList& regionNames = prePacing.regionNames();
    List<scalarField> regionStates(regionNames.size());
    scalarField regionTa(regionNames.size(), 0.0);
    labelList beats(regionNames.size(), 0);

    forAll(regionNames, regionI)
    {
        if
        (
            !prePacing.paced(regionI)
         || Pstream::myProcNo() != myocardiumPrePacing::owner(regionI)
        )
        {
            continue;
        }

        const prePacingIO::PrePacingConfig& cfg = prePacing.config(regionI);
        autoPtr<ionicModel> cellPtr = prePacing.newSingleCell(regionI, true);

        autoPtr<activeTensionModel> tensionPtr =
            activeTensionModel::New(electroMechanicalProperties(), 1);
        tensionPtr->setElectromechanicalSignalProvider(cellPtr());
        tensionPtr->validateProvider();

        if (!tensionPtr->states())
        {
            FatalErrorInFunction
                << "Active tension model '" << tensionPtr->type()
                << "' does not expose its states; prePacing cannot pace it."
                << exit(FatalError);
        }

        isometricTension coupling(tensionPtr());

        beats[regionI] =
            cellPtr->prePaceToConvergence
            (
                prePacing.singleCellDt(regionI),
                cfg.tolerance,
                cfg.minBeats,
                cfg.maxBeats,
                cfg.beatComparisonInterval,
                &coupling
            );

        regionStates[regionI] = coupling.state();
        regionTa[regionI] = coupling.Ta();
    }

    // Share each region's result from the processor that paced it.
    struct takeComputed
    {
        void operator()(scalarField& x, const scalarField& y) const
        {
            if (x.empty())
            {
                x = y;
            }
        }
    };
    Pstream::listCombineReduce(regionStates, takeComputed());
    Pstream::listCombineReduce(beats, maxEqOp<label>());
    Pstream::listCombineReduce(regionTa, plusEqOp<scalar>());

    List<scalarField> states(*tissueStatesPtr);
    labelList regionCellCount(regionNames.size(), 0);

    forAll(cellRegion, cellI)
    {
        const label regionI = cellRegion[cellI];
        if (prePacing.paced(regionI))
        {
            states[cellI] = regionStates[regionI];
            ++regionCellCount[regionI];
        }
    }

    if (!activeTensionModel_->setStates(states))
    {
        FatalErrorInFunction
            << "Active tension model '" << activeTensionModel_->type()
            << "' does not accept seeded states; prePacing cannot seed it."
            << exit(FatalError);
    }

    forAll(regionNames, regionI)
    {
        if (prePacing.paced(regionI))
        {
            Info<< "prePacing: active tension of region '"
                << regionNames[regionI] << "' converged after "
                << beats[regionI] << " beats; seeded "
                << regionCellCount[regionI] << " cells, diastolic Ta = "
                << regionTa[regionI] << " (model units)." << endl;
        }
    }
}


void sequentialElectroMechanical::updateLambda()
{
    const fvMesh& solidMesh = solid().mesh();

    const bool hasD  = solidMesh.foundObject<volVectorField>("D");
    const bool hasF0 = solidMesh.foundObject<volVectorField>("f0");

    if (!hasD || !hasF0)
    {
        if (activeTensionRequirements_.needsLambda)
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field "
                << (!hasD ? "D" : "f0")
                << " disappeared from the solid objectRegistry at t="
                << runTime().value() << "."
                << abort(FatalError);
        }
        lambdaField_ = 1.0;
        return;
    }

    const volVectorField& D  = solidMesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = solidMesh.lookupObject<volVectorField>("f0");

    // Deformation gradient F = I + grad(D)^T (total Lagrangian convention,
    // matching solids4foam mechanicalLaw). The fibre stretch follows from
    // lambda^2 = f0 & C & f0 = (F & f0) & (F & f0), i.e. lambda = mag(F & f0).
    const volTensorField gradD(fvc::grad(D));

    forAll(lambdaField_, cellI)
    {
        const tensor F(I + gradD[cellI].T());
        lambdaField_[cellI] = mag(F & f0[cellI]);
    }
}


bool sequentialElectroMechanical::evolve()
{
    Info<< "Evolving " << type() << endl;

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->preSolve
        (
            const_cast<volScalarField&>(electro().Vm()),
            solid().D()
        );
    }

    electro().evolve();

    // Update the fibre stretch from the (lagged) solid deformation before
    // evaluating the active tension.
    updateLambda();

    const scalar t  = runTime().value();
    const scalar dt = runTime().deltaT().value();

    scalarField& TaI = Ta_.primitiveFieldRef();

    activeTensionModel_->calculateTension(t, dt, lambdaField_, TaI);

    if (TaScale_ != 1.0)  // skip no-op multiply; 1.0 is exactly representable
    {
        TaI *= TaScale_;
    }

    Ta_.correctBoundaryConditions();

    solid().evolve();
    solid().updateTotalFields();

    if
    (
        verificationModelPtr_.valid()
     && verificationModelPtr_->shouldPostProcess(electro().Vm(), solid().D())
    )
    {
        verificationModelPtr_->postProcess(electro().Vm(), solid().D(), Ta_);
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// ************************************************************************* //
