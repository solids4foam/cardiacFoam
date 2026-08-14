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
#include "fixedValueFvPatchFields.H"
#include "conductivityFieldIO.H"
#include "eikonalVerificationModel.H"
#include "zeroGradientFvPatchFields.H"

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


// NOTE: a free function here used to re-derive "is verification enabled?"
// from type + enabled, in parallel with eikonalVerificationModel::New, which
// answers the same question by returning nullptr. The boundary types keyed
// off that duplicate rather than off the verifier that actually got built.
// The single source of truth is now verificationModelPtr_ itself.


// activationTime's boundary conditions come from 0/activationTime like any
// other field's. They used to be chosen here from the verification settings
// and handed to the field constructor, which hid them from the case: nothing
// in the case directory said whether the run used zeroGradient or fixedValue.
//
// Manufactured-solution runs need Dirichlet boundaries -- applyConstraints
// writes the exact solution onto each patch face, and a zeroGradient patch
// would extrapolate those values away from the interior on the next
// evaluate(), silently invalidating the verification while still reporting
// error norms. checkVerificationBoundaryTypes below enforces that instead.
void checkVerificationBoundaryTypes(const volScalarField& activationTime)
{
    const volScalarField::Boundary& boundary = activationTime.boundaryField();

    forAll(boundary, patchI)
    {
        const fvPatchScalarField& patchField = boundary[patchI];

        if (patchField.empty() || patchField.type() == "empty")
        {
            continue;
        }

        if (!patchField.fixesValue())
        {
            FatalErrorInFunction
                << "Manufactured-solution verification is active, but patch '"
                << patchField.patch().name() << "' of activationTime is of "
                << "type '" << patchField.type() << "', which does not fix "
                << "its value." << nl
                << "The verifier imposes the exact solution on the boundary; "
                << "a non-Dirichlet patch discards it and the reported error "
                << "norms become meaningless." << nl
                << "Set this patch to fixedValue in 0/activationTime."
                << exit(FatalError);
        }
    }
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


tmp<volTensorField> eikonalMyocardiumDomain::initialiseConductivity() const
{
    return readConductivityField
    (
        mesh(),
        supportMesh_,
        meshSubsetPtr_.valid() ? &meshSubsetPtr_() : nullptr,
        electroProperties_,
        conductivityFieldSpec
        {
            "Conductivity",
            "conductivity",
            "conductivity"
        }
    );
}


eikonalMyocardiumDomain::eikonalMyocardiumDomain
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
            "activationTime",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_)
    ),
    Vm_
    (
        IOobject
        (
            "Vm",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            max
            (
                gradActivationTime_ & w_,
                dimensionedScalar("zero", dimTime, 0.0)
            )
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
    ),
    verificationModelPtr_()
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


    verificationModelPtr_ =
        eikonalVerificationModel::New
        (
            electroProperties_,
            mesh(),
            conductivity_,
            chi_,
            Cm_,
            c0_,
            eikonalAdvectionDiffusionApproach_
        );

    if (verificationModelPtr_)
    {
        checkVerificationBoundaryTypes(activationTime_);
    }
}


eikonalMyocardiumDomain::~eikonalMyocardiumDomain() = default;


void eikonalMyocardiumDomain::preInitialiseFromSeeds
(
    const labelList& constrainedCells,
    const scalarField& constrainedValues
)
{
    scalarField& T = activationTime_.primitiveFieldRef();

    // Upper-bound CV from the fastest propagation direction of M
    scalar maxMdiag = 0;
    forAll(M_, cellI)
    {
        maxMdiag = max(maxMdiag, max(M_[cellI].xx(), max(M_[cellI].yy(), M_[cellI].zz())));
    }
    reduce(maxMdiag, maxOp<scalar>());
    const scalar CV_est = c0_.value() * Foam::sqrt(maxMdiag + SMALL);

    // Seed cells keep their constrained values; all others go to GREAT
    forAll(T, cellI) { if (T[cellI] < 0) T[cellI] = GREAT; }
    activationTime_.correctBoundaryConditions();

    // Bellman-Ford relay: propagate minimum arrival time across face connectivity.
    // Each pass reduces T for cells reachable from seeds; parallel-safe because
    // correctBoundaryConditions syncs processor patches after every pass.
    const vectorField& cc = mesh().cellCentres();
    const label nIntFaces  = mesh().nInternalFaces();
    const labelList& own   = mesh().faceOwner();
    const labelList& nei   = mesh().faceNeighbour();

    bool changed = true;
    while (changed)
    {
        changed = false;
        for (label faceI = 0; faceI < nIntFaces; ++faceI)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            const scalar dt = mag(cc[o] - cc[n]) / CV_est;

            if (T[o] < GREAT/2 && T[o] + dt < T[n]) { T[n] = T[o] + dt; changed = true; }
            if (T[n] < GREAT/2 && T[n] + dt < T[o]) { T[o] = T[n] + dt; changed = true; }
        }

        forAll(activationTime_.boundaryField(), patchI)
        {
            const fvPatchScalarField& pT = activationTime_.boundaryField()[patchI];
            if (pT.coupled())
            {
                const coupledFvPatchScalarField& cpT = refCast<const coupledFvPatchScalarField>(pT);
                tmp<scalarField> tneiT = cpT.patchNeighbourField();
                const scalarField& nT = tneiT();

                const labelUList& faceCells = pT.patch().faceCells();
                tmp<vectorField> tdelta = pT.patch().delta();
                const vectorField& delta = tdelta();

                forAll(faceCells, faceI)
                {
                    const label cellI = faceCells[faceI];
                    const scalar dt = mag(delta[faceI]) / CV_est;

                    if (nT[faceI] < GREAT/2 && nT[faceI] + dt < T[cellI])
                    {
                        T[cellI] = nT[faceI] + dt;
                        changed = true;
                    }
                }
            }
        }

        activationTime_.correctBoundaryConditions();
        reduce(changed, orOp<bool>());
    }

    // Re-enforce exact seed values (BF may have relaxed them from a closer seed)
    forAll(constrainedCells, i) { T[constrainedCells[i]] = constrainedValues[i]; }

    // Cells with no path to any seed (disconnected regions) → restore sentinel
    forAll(T, cellI) { if (T[cellI] >= GREAT/2) T[cellI] = -1; }

    activationTime_.correctBoundaryConditions();
}


void eikonalMyocardiumDomain::advance
(
    scalar t0,
    scalar dt
)
{
    advance(t0, dt, nullptr);
}


void eikonalMyocardiumDomain::advance
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    (void)t0;
    (void)dt;

    scalarField& activationValues = activationTime_.primitiveFieldRef();
    forAll(stimulusCellIDs_, i)
    {
        activationValues[stimulusCellIDs_[i]] = 0.0;
    }
    activationTime_.correctBoundaryConditions();

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->applyConstraints(activationTime_);
    }

    labelList constrainedCells;
    scalarField constrainedValues;
    collectActivationConstraints
    (
        activationTime_,
        constrainedCells,
        constrainedValues
    );

    preInitialiseFromSeeds(constrainedCells, constrainedValues);

    const dimensionedScalar one("one", dimless, 1.0);
    const dimensionedScalar smallG("smallG", dimTime, SMALL);

    tmp<volScalarField> tSmms;
    if (verificationModelPtr_.valid() && verificationModelPtr_->enabled())
    {
        tSmms = verificationModelPtr_->sourceTerm();
    }
    else
    {
        tSmms.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Smms",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("Smms", dimless, 0.0)
            )
        );
    }
    const volScalarField& Smms = tSmms();

    auto solveActivationEqn = [&]()
    {
        gradActivationTime_ = fvc::grad(activationTime_);
        w_ = M_ & gradActivationTime_;
        // Guard against negative dot product: M is positive-definite so
        // dot(grad, M*grad) >= 0 analytically, but NaN/Inf from a stalled
        // inner solve or degenerate mesh cells can violate this numerically.
        G_ = sqrt(max(gradActivationTime_ & w_, dimensionedScalar("zero", dimTime, 0.0)) + smallG);

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
              + Smms
            );

            // Under-relax the nonlinear outer loop: G, phiU and divPhiU are
            // lagged (recomputed from psi each outer iteration), so without
            // relaxation the deferred-correction fixed point is not contractive
            // and the outer residual drifts upwards.  relax() also boosts
            // diagonal dominance, which stabilises the unpreconditioned solve.
            activationEqn.relax();
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
              + Smms
            );

            activationEqn.relax();
            activationEqn.setValues(constrainedCells, constrainedValues);
            activationEqn.solve();
        }
    };

    if (pimplePtr)
    {
        while (pimplePtr->loop())
        {
            solveActivationEqn();
        }
    }
    else
    {
        solveActivationEqn();
    }
}


scalar eikonalMyocardiumDomain::suggestExplicitDeltaT(scalar maxCo) const
{
    (void)maxCo;
    return 1.0;
}


bool eikonalMyocardiumDomain::applyModelTimeControls(Time& runTime)
{
    InfoInFunction << "Setting deltaT and endTime to 1.0" << endl;
    runTime.setDeltaT(1.0);
    runTime.setEndTime(runTime.deltaT());
    return true;
}


bool eikonalMyocardiumDomain::shouldPostProcess() const
{
    return
        verificationModelPtr_.valid()
     && verificationModelPtr_->shouldPostProcess();
}


void eikonalMyocardiumDomain::postProcess()
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->postProcess(activationTime_);
    }
}


void eikonalMyocardiumDomain::write()
{
    if (meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh())
    {
        const labelUList& cellMap = meshSubsetPtr_->cellMap();
        writeMappedCellField(activationTime_, supportMesh_, cellMap);
        return;
    }

    activationTime_.write();
}

} // End namespace Foam

// ************************************************************************* //
