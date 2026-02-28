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

#include "monoDomainElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(monoDomainElectro, 0);
addToRunTimeSelectionTable(electroModel, monoDomainElectro, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<volTensorField> monoDomainElectro::initialiseConductivity() const
{
    // Prepare the result
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "conductivity",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor
            (
                "zero",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                tensor::zero
            )
        )
    );
    volTensorField& result = tresult.ref();

    // Check if we actually read it, otherwise assign constant value
    if (!result.headerOk())
    {
        Info<< "\nconductivity not found on disk, using conductivity from "
            << cardiacProperties_.name() << nl << endl;

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "conductivity",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    cardiacProperties_
                ) & tensor(I)
            );

        if (result.size() > 0)
        {
            Info<< "Conductivity tensor (cell 0): " << result[0] << nl;
        }
    }
    else
    {
        Info<< "conductivity field read from " << runTime().timeName() << nl
            << endl;
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

monoDomainElectro::monoDomainElectro
(
    Time& runTime,
    const word& region
)
:
    electroModel(typeName, runTime, region),
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
        dimensionedScalar("Vm", dimVoltage, -0.084),
        "zeroGradient"
    ),
    Iion_
    (
        IOobject
        (
            "ionicCurrent",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVoltage/ dimTime, 0.0),
        "zeroGradient"
    ),
    externalStimulusCurrent_
    (
        IOobject
        (
            "externalStimulusCurrent",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    externalStimulus_(),
    cardiacProperties_(electroProperties().subDict("cardiacProperties")),
    chi_("chi", dimArea/dimVolume, cardiacProperties_),
    Cm_("Cm",dimCurrent*dimTime/(dimVoltage*dimArea), cardiacProperties_),
    conductivity_(initialiseConductivity()),
    ionicModelPtr_
    (
        electroProperties().subDict("ionicModel"),
        mesh().nCells(),
        runTime.deltaTValue()
    )
{
    // Initialise the external stimulii
    const dictionary& stimulusProtocol =
        electroProperties().subDict("stimulusProtocol");

    // Read external stimulus from one or more bounding boxes
    List<boundBox> stimulusBoxes;
    List<scalar> stimulusStartTimes;

    if
    (
        stimulusProtocol.found("stimulusLocationMinList")
     && stimulusProtocol.found("stimulusLocationMaxList")
    )
    {
        const List<point> mins
        (
            stimulusProtocol.lookup("stimulusLocationMinList")
        );
        const List<point> maxs
        (
            stimulusProtocol.lookup("stimulusLocationMaxList")
        );

        if (mins.size() != maxs.size())
        {
            FatalErrorInFunction
                << "stimulusLocationMinList and stimulusLocationMaxList must "
                << "have the same size."
                << abort(FatalError);
        }

        stimulusBoxes.setSize(mins.size());
        forAll(mins, i)
        {
            stimulusBoxes[i] = boundBox(mins[i], maxs[i]);
        }
    }
    else
    {
        stimulusBoxes.setSize(1);
        stimulusBoxes[0] = boundBox
        (
            point(stimulusProtocol.lookup("stimulusLocationMin")),
            point(stimulusProtocol.lookup("stimulusLocationMax"))
        );
    }

    if (stimulusProtocol.found("stimulusStartTimeList"))
    {
        stimulusProtocol.lookup("stimulusStartTimeList")
            >> stimulusStartTimes;

        if (stimulusStartTimes.size() != stimulusBoxes.size())
        {
            FatalErrorInFunction
                << "stimulusStartTimeList must have the same size as "
                << "the stimulus box list."
                << abort(FatalError);
        }
    }
    else
    {
        stimulusStartTimes.setSize(stimulusBoxes.size());
        const scalar startTime =
            stimulusProtocol.lookupOrDefault<scalar>("stimulusStartTime", 0.0);
        forAll(stimulusStartTimes, i)
        {
            stimulusStartTimes[i] = startTime;
        }
    }

    List<labelHashSet> stimCellSets(stimulusBoxes.size());
    labelHashSet stimCellSet;

    forAll(mesh.C(), cellI)
    {
        const point& c = mesh.C()[cellI];
        forAll(stimulusBoxes, bI)
        {
            if (stimulusBoxes[bI].contains(c))
            {
                stimCellSets[bI].insert(cellI);
                stimCellSet.insert(cellI);
                break;
            }
        }
    }

    const dimensionedScalar stimulusDuration
    (
        "stimulusDuration", dimTime, stimulusProtocol
    );

    const dimensionedScalar stimulusIntensity
    (
        "stimulusIntensity", dimCurrent/dimVolume, stimulusProtocol
    );

    // Set the number of stimulii
    externalStimulus_.setSize(mins.size());

    // Initialise each stimulii
    forAll(stimulusBoxes, bI)
    {
        externalStimulus_[bI].cellIDs_ = stimCellSets[bI].toc();
        externalStimulus_[bI].startTime_ = stimulusStartTimes[bI];
        externalStimulus_[bI].duration_ = stimulusDuration;
        externalStimulus_[bI].stimulusIntensity_ = stimulusIntensity;
    }

    // Only print stimulus info if this is a 3D mesh
    if (mesh.nGeometricD() == 3)
    {
        Info<< "---------------------------------------------" << nl
            << "3D mesh detected — stimulus parameters:"       << nl
            << "Stimulus boxes: " << stimulusBoxes.size()      << nl
            << "Number of cells in stimulus region: "
            << stimCellSet.size()                              << nl
            << "Stimulus duration: "
            << stimulusDuration.value() << " "
            << stimulusDuration.dimensions()                   << nl
            << "Stimulus start time(s): " << stimulusStartTimes << nl
            << "Stimulus intensity: "
            << stimulusIntensity.value() << " "
            << stimulusIntensity.dimensions()                  << nl
            << "---------------------------------------------" << endl;
    }

    // Write parameters
    Info<< "Surface-to-volume ratio chi = " << chi_ << nl
        << "Membrane capacitance Cm = " << Cm_ << nl
        << "M tensor field:" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool monoDomainElectro::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    /*
    // Define convenient dimensioned scalars
    const dimensionedScalar one("one", dimless, 1.0);
    const dimensionedScalar smallG("smallG", dimTime, SMALL);

    while (pimple().loop())
    {
        // Update the gradient
        gradPsi_ = fvc::grad(psi_);

        // Update w term
        w_ = M_ & gradPsi_;

        // Update the nonlinear term
        G_ = sqrt((gradPsi_ & w_) + smallG);

        // Solve eikonal equation
        if (eikonalAdvectionDiffusionApproach_)
        {
            // Update a term
            a_ = w_/G_;

            // Update the u term
            u_ = c0_*a_;

            // Update the phiU term
            phiU_ = (fvc::interpolate(u_) & mesh().Sf());

            // Update div(phiU)
            divPhiU_ = fvc::div(phiU_);

            // Construct the eikonal–diffusion equation (Eq. 6.14 in Quarteroni
            // et al.)
            //
            // The equation is written in a stabilised Picard (defect-
            // correction) form.
            // The implicit advection operator on the LHS is derived from a
            // linearisation of the nonlinear eikonal term c0*G, but is used
            // here purely to improve convergence. The corresponding explicit
            // advection term on the RHS ensures that, at convergence of the
            // outer iterations, the original eikonal–diffusion equation is
            // recovered exactly.
            //
            // The SuSp formulation is used to avoid diagonal weakening arising from
            // div(phiU), improving robustness of the linear solve.
            fvScalarMatrix psiEqn
            (
              - fvm::laplacian(M_, psi_)
              + fvm::div(phiU_, psi_) + fvm::SuSp(-divPhiU_, psi_)
             == one
              + fvc::div(phiU_, psi_) - divPhiU_*psi_
              - c0_*G_
            );

            // Enforce psi to be zero for the stimulus cells
            psiEqn.setValues(stimulusCellIDs_, 0.0);

            // Solve the system for psi
            psiEqn.solve("asymmetric_" + psi_.name());
        }
        else
        {
            // Construct the eikonal-diffusion equation, where the nonlinear
            // term c0*G is treated entirely as an explicit deferred correction
            fvScalarMatrix psiEqn
            (
              - fvm::laplacian(M_, psi_)
              + c0_*G_
             == one
            );

            // Optional equation relaxation
            psiEqn.relax();

            // Enforce psi to be zero for the stimulus cells
            psiEqn.setValues(stimulusCellIDs_, 0.0);

            // Solve the system for psi
            psiEqn.solve();
        }
    }
        */
        
    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
