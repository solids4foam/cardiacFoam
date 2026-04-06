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

#include "eikonalSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// OverrideTypeName("eikonalSolver") in the header declares typeName_() == "eikonalSolver".
// defineTypeNameWithName registers the static member accordingly; the plain
// defineTypeNameAndDebug(EikonalSolver, 0) would use #EikonalSolver and overwrite it.
defineTypeNameWithName(EikonalSolver, "eikonalSolver");
defineDebugSwitch(EikonalSolver, 0);
addToRunTimeSelectionTable(electroModel, EikonalSolver, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<volTensorField> EikonalSolver::initialiseConductivity() const
{
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

    if (!result.headerOk())
    {
        Info<< "\nconductivity not found on disk, using conductivity from "
            << electroProperties().name() << nl << endl;

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "conductivity",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    electroProperties()
                )
              & tensor(I)
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

EikonalSolver::EikonalSolver(Time& runTime, const word& region)
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
    gradPsi_(fvc::grad(psi_)),
    stimulusCellIDs_(0),
    chi_("chi", dimArea/dimVolume, electroProperties()),
    Cm_("cm", dimCurrent*dimTime/(dimVoltage*dimArea), electroProperties()),
    conductivity_(initialiseConductivity()),
    M_("M", conductivity_/(chi_*Cm_)),
    w_("w", M_ & gradPsi_),
    G_
    (
        "G",
        sqrt((gradPsi_ & w_) + dimensionedScalar("smallG", dimTime, SMALL))
    ),
    c0_("c0", electroProperties()),
    a_("a", w_/G_),
    u_("u", c0_*a_),
    phiU_("phiU", (fvc::interpolate(u_) & mesh().Sf())),
    divPhiU_("divPhiU", fvc::div(phiU_)),
    eikonalAdvectionDiffusionApproach_
    (
        electroProperties().lookup("eikonalAdvectionDiffusionApproach")
    )
{
    const boundBox bb
    (
        point(electroProperties().lookup("stimulusLocationMin")),
        point(electroProperties().lookup("stimulusLocationMax"))
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

    stimulusCellIDs_ = stimCellSet.toc();

    Info<< "Surface-to-volume ratio chi = " << chi_ << nl
        << "Membrane capacitance Cm = " << Cm_ << nl
        << "M tensor field:" << nl
        << "    max(M) = " << gMax(M_) << nl
        << "    min(M) = " << gMin(M_) << nl
        << "    average(M) = " << gAverage(M_) << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool EikonalSolver::evolve()
{
    Info<< "Evolving electro model: " << this->type() << endl;

    const dimensionedScalar one("one", dimless, 1.0);
    const dimensionedScalar smallG("smallG", dimTime, SMALL);

    while (pimple().loop())
    {
        gradPsi_ = fvc::grad(psi_);
        w_ = M_ & gradPsi_;
        G_ = sqrt((gradPsi_ & w_) + smallG);

        if (eikonalAdvectionDiffusionApproach_)
        {
            a_ = w_/G_;
            u_ = c0_*a_;
            phiU_ = (fvc::interpolate(u_) & mesh().Sf());
            divPhiU_ = fvc::div(phiU_);

            fvScalarMatrix psiEqn
            (
               -fvm::laplacian(M_, psi_)
              + fvm::div(phiU_, psi_)
              + fvm::SuSp(-divPhiU_, psi_)
             == one + fvc::div(phiU_, psi_) - divPhiU_*psi_ - c0_*G_
            );

            psiEqn.setValues(stimulusCellIDs_, 0.0);
            psiEqn.solve("asymmetric_" + psi_.name());
        }
        else
        {
            fvScalarMatrix psiEqn
            (
               -fvm::laplacian(M_, psi_)
              + c0_*G_
             == one
            );

            psiEqn.relax();
            psiEqn.setValues(stimulusCellIDs_, 0.0);
            psiEqn.solve();
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
