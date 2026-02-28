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

#include "eikonalDiffusionElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eikonalDiffusionElectro, 0);
addToRunTimeSelectionTable(electroModel, eikonalDiffusionElectro, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<volTensorField> eikonalDiffusionElectro::initialiseConductivity() const
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

eikonalDiffusionElectro::eikonalDiffusionElectro
(
    Time& runTime,
    const word& region
)
:
    electroModel(typeName, runTime, region),
    psi_
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("psi", dimTime, 0.0),
        "zeroGradient"
    ),
    gradPsi_(fvc::grad(psi_)),
    stimulusCellIDs_(0),
    cardiacProperties_(electroProperties().subDict("cardiacProperties")),
    chi_("chi", dimArea/dimVolume, cardiacProperties_),
    Cm_("Cm",dimCurrent*dimTime/(dimVoltage*dimArea), cardiacProperties_),
    conductivity_(initialiseConductivity()),
    M_("M", conductivity_/(chi_*Cm_)),
    w_("w", M_ & gradPsi_),
    G_("G", sqrt((gradPsi_ & w_) + dimensionedScalar("smallG", dimTime, SMALL))),
    c0_("c0", cardiacProperties_),
    a_("a", w_/G_),
    u_("u", c0_*a_),
    phiU_("phiU", (fvc::interpolate(u_) & mesh().Sf())),
    divPhiU_("divPhiU", fvc::div(phiU_)),
    eikonalAdvectionDiffusionApproach_
    (
        electroProperties().lookup("eikonalAdvectionDiffusionApproach")
    )
{
    // Initialise stimulusCellIDs
    Info<< "Reading stimulus protocol\n" << endl;
    const dictionary& stimulusProtocol
    (
        electroProperties().subDict("stimulusProtocol")
    );

    // Read external stimulus from a bounding box
    // Find cells in the stimulus volume
    const boundBox bb
    (
        point(stimulusProtocol.lookup("stimulusLocationMin")),
        point(stimulusProtocol.lookup("stimulusLocationMax"))
    );

    labelHashSet stimCellSet;

    const fvMesh& mesh = this->mesh();
    forAll(mesh.C(), cellI)
    {
        if (bb.contains(mesh.C()[cellI]))
        {
            stimCellSet.insert(cellI);
        }
    }

    // Record the simulus cells
    stimulusCellIDs_ = stimCellSet.toc();

    // Write parameters
    Info<< "Surface-to-volume ratio chi = " << chi_ << nl
        << "Membrane capacitance Cm = " << Cm_ << nl
        << "M tensor field:" << nl
        << "    max(M) = " << gMax(M_) << nl
        << "    min(M) = " << gMin(M_) << nl
        << "    average(M) = " << gAverage(M_) << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool eikonalDiffusionElectro::evolve()
{
    Info<< "Evolving electro model: " << this->type() << endl;

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

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
