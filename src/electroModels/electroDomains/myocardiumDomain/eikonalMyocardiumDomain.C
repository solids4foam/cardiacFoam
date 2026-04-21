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

#include "eikonalMyocardiumDomain.H"
#include "error.H"
#include "DynamicList.H"
#include "Switch.H"

namespace Foam
{

namespace
{

autoPtr<fvMeshSubset> createMyocardiumMeshSubset
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
)
{
    if (!electroProperties.found("cellZone"))
    {
        return autoPtr<fvMeshSubset>(nullptr);
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = supportMesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << supportMesh.name() << "'."
            << exit(FatalError);
    }

    autoPtr<fvMeshSubset> subsetPtr(new fvMeshSubset(supportMesh));
    subsetPtr->setCellSubset(supportMesh.cellZones()[zoneId]);
    return subsetPtr;
}


const fvMesh& resolveMyocardiumMesh
(
    const fvMesh& supportMesh,
    const autoPtr<fvMeshSubset>& subsetPtr
)
{
    return
    (
        subsetPtr.valid() && subsetPtr->hasSubMesh()
      ? subsetPtr->subMesh()
      : supportMesh
    );
}


void writeMappedCellField
(
    const volScalarField& subField,
    const fvMesh& supportMesh,
    const labelUList& cellMap
)
{
    volScalarField fullField
    (
        IOobject
        (
            subField.name(),
            supportMesh.time().timeName(),
            supportMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        supportMesh,
        dimensionedScalar("zero", subField.dimensions(), 0.0),
        "zeroGradient"
    );

    scalarField& fullValues = fullField.primitiveFieldRef();
    const scalarField& subValues = subField.primitiveField();

    forAll(cellMap, subCelli)
    {
        fullValues[cellMap[subCelli]] = subValues[subCelli];
    }

    fullField.write();
}


void collectActivationConstraints
(
    const volScalarField& activationTime,
    labelList& constrainedCells,
    scalarField& constrainedValues
)
{
    DynamicList<label> cells;
    DynamicList<scalar> values;

    const scalarField& activationValues = activationTime.primitiveField();

    forAll(activationValues, cellI)
    {
        if (activationValues[cellI] >= 0.0)
        {
            cells.append(cellI);
            values.append(activationValues[cellI]);
        }
    }

    constrainedCells = cells;
    constrainedValues = values;
}

} // End anonymous namespace


tmp<volTensorField> EikonalMyocardiumDomain::initialiseConductivity() const
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "conductivity",
                mesh().time().timeName(),
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
        if (electroProperties_.lookupOrDefault<Switch>("reportSetup", false))
        {
            Info<< "\nconductivity not found on disk, using conductivity from "
                << electroProperties_.name()
                << nl << endl;
        }

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "conductivity",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    electroProperties_
                )
              & tensor(I)
            );
    }

    return tresult;
}


EikonalMyocardiumDomain::EikonalMyocardiumDomain
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
)
:
    meshSubsetPtr_(createMyocardiumMeshSubset(supportMesh, electroProperties)),
    supportMesh_(supportMesh),
    activationTime_
    (
        IOobject
        (
            "psi",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("psi", dimTime, -1.0),
        "zeroGradient"
    ),
    Vm_
    (
        IOobject
        (
            "Vm",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("Vm", dimVoltage, -80.0),
        "zeroGradient"
    ),
    sourceField_
    (
        IOobject
        (
            "externalStimulusCurrent",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    gradActivationTime_(fvc::grad(activationTime_)),
    stimulusCellIDs_(0),
    electroProperties_(electroProperties),
    chi_("chi", dimArea/dimVolume, electroProperties),
    Cm_("cm", dimCurrent*dimTime/(dimVoltage*dimArea), electroProperties),
    conductivity_(initialiseConductivity()),
    M_("M", conductivity_/(chi_*Cm_)),
    w_("w", M_ & gradActivationTime_),
    G_
    (
        "G",
        sqrt
        (
            (gradActivationTime_ & w_)
          + dimensionedScalar("smallG", dimTime, SMALL)
        )
    ),
    c0_("c0", electroProperties),
    a_("a", w_/G_),
    u_("u", c0_*a_),
    phiU_("phiU", (fvc::interpolate(u_) & mesh().Sf())),
    divPhiU_("divPhiU", fvc::div(phiU_)),
    eikonalAdvectionDiffusionApproach_
    (
        electroProperties.lookup("eikonalAdvectionDiffusionApproach")
    )
{
    const boundBox bb
    (
        point(electroProperties.lookup("stimulusLocationMin")),
        point(electroProperties.lookup("stimulusLocationMax"))
    );

    labelHashSet stimCellSet;
    const fvMesh& myocardiumMesh = mesh();

    forAll(myocardiumMesh.C(), cellI)
    {
        if (bb.contains(myocardiumMesh.C()[cellI]))
        {
            stimCellSet.insert(cellI);
        }
    }

    stimulusCellIDs_ = stimCellSet.toc();

    if
    (
        electroProperties.lookupOrDefault<Switch>("reportSetup", false)
     && meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh()
    )
    {
        Info<< "Constructed EikonalMyocardiumDomain on submesh '"
            << mesh().name() << "' from cellZone '"
            << electroProperties.lookupOrDefault<word>("cellZone", word::null)
            << "'." << nl << endl;
    }
}


void EikonalMyocardiumDomain::advance
(
    scalar t0,
    scalar dt
)
{
    advance(t0, dt, nullptr);
}


void EikonalMyocardiumDomain::advance
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    (void)t0;
    (void)dt;
    (void)pimplePtr;

    scalarField& activationValues = activationTime_.primitiveFieldRef();
    forAll(stimulusCellIDs_, i)
    {
        activationValues[stimulusCellIDs_[i]] = 0.0;
    }
    activationTime_.correctBoundaryConditions();

    labelList constrainedCells;
    scalarField constrainedValues;
    collectActivationConstraints
    (
        activationTime_,
        constrainedCells,
        constrainedValues
    );

    const dimensionedScalar one("one", dimless, 1.0);
    const dimensionedScalar smallG("smallG", dimTime, SMALL);

    gradActivationTime_ = fvc::grad(activationTime_);
    w_ = M_ & gradActivationTime_;
    G_ = sqrt((gradActivationTime_ & w_) + smallG);

    if (eikonalAdvectionDiffusionApproach_)
    {
        a_ = w_/G_;
        u_ = c0_*a_;
        phiU_ = (fvc::interpolate(u_) & mesh().Sf());
        divPhiU_ = fvc::div(phiU_);

        fvScalarMatrix activationEqn
        (
           -fvm::laplacian(M_, activationTime_)
          + fvm::div(phiU_, activationTime_)
          + fvm::SuSp(-divPhiU_, activationTime_)
         == one
          + fvc::div(phiU_, activationTime_)
          - divPhiU_*activationTime_
          - c0_*G_
        );

        activationEqn.setValues(constrainedCells, constrainedValues);
        activationEqn.solve("asymmetric_" + activationTime_.name());
    }
    else
    {
        fvScalarMatrix activationEqn
        (
           -fvm::laplacian(M_, activationTime_)
          + c0_*G_
         == one
        );

        activationEqn.relax();
        activationEqn.setValues(constrainedCells, constrainedValues);
        activationEqn.solve();
    }
}


scalar EikonalMyocardiumDomain::suggestExplicitDeltaT(scalar maxCo) const
{
    (void)maxCo;
    return 1.0;
}


bool EikonalMyocardiumDomain::applyModelTimeControls(Time& runTime) const
{
    InfoInFunction << "Setting deltaT and endTime to 1.0" << endl;
    runTime.setDeltaT(1.0);
    runTime.setEndTime(runTime.deltaT());
    return true;
}


void EikonalMyocardiumDomain::write()
{
    if (meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh())
    {
        const labelUList& cellMap = meshSubsetPtr_->cellMap();

        writeMappedCellField(activationTime_, supportMesh_, cellMap);
        writeMappedCellField(Vm_, supportMesh_, cellMap);
        writeMappedCellField(sourceField_, supportMesh_, cellMap);
        return;
    }

    activationTime_.write();
    Vm_.write();
    sourceField_.write();
}

} // End namespace Foam

// ************************************************************************* //
