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
#include "insulatedFaceConductivity.H"
#include "sealedWallPolicy.H"

#include <cmath>
#include <limits>

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
    subsetPtr->setCellSubset
    (
        supportMesh.cellZones()[zoneId],
        myocardiumExposedFacesPatch(supportMesh, electroProperties)
    );
#if OPENFOAM < 2406
    // OpenFOAM before 2406 copies the base fvSchemes and fvSolution onto the
    // sub-mesh without parsing them, so read them for the scheme and solver
    // lookups to resolve.
    fvMesh& subMesh = subsetPtr->subMesh();
    subMesh.fvSchemes::readOpt(IOobject::MUST_READ);
    subMesh.fvSchemes::read();
    subMesh.fvSolution::readOpt(IOobject::MUST_READ);
    subMesh.fvSolution::read();
#endif
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


// activationTime's boundary conditions come from 0/activationTime like any
// other field's.
//
// Manufactured-solution runs need Dirichlet boundaries -- applyConstraints
// writes the exact solution onto each patch face, and a zeroGradient patch
// would extrapolate those values away from the interior on the next
// evaluate(), silently invalidating the verification while still reporting
// error norms. checkVerificationBoundaryTypes below enforces that.
void checkVerificationBoundaryTypes
(
    const volScalarField& activationTime,
    const eikonalVerificationModel& verifier
)
{
    const volScalarField::Boundary& boundary = activationTime.boundaryField();

    forAll(boundary, patchI)
    {
        const fvPatchScalarField& patchField = boundary[patchI];

        if
        (
            patchField.empty()
         || patchField.type() == "empty"
         || patchField.patch().coupled()
         || !verifier.imposesExactValue(patchField.patch())
        )
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


// Fills one field of the metric G = sqrt(Q), Q = g & (M & g), and of the
// direction a = (M & g)/G, given w = M & g.
//
// Q is non-negative for positive-definite M. Evaluated in floating point
// (the matrix-vector product w, then the inner product g & w) it differs from
// the exact value by at most gamma_6 (|g| & (|M| & |g|)), with
// gamma_k = k u/(1 - k u) and u the unit roundoff (Higham, Accuracy and
// Stability of Numerical Algorithms, 2nd ed., Sec. 3.1 and 3.5). A negative Q
// within that bound is rounding of a zero or vanishing Q and is set to zero;
// a negative Q beyond it means M is not positive definite, and a non-finite Q
// means g or M is not finite: both are fatal.
//
// a = w/G where Q > 0 and a = 0 where Q = 0, so a & g = G everywhere. For
// symmetric positive-definite M, |a| = |M^1/2 (M^1/2 g)|/|M^1/2 g|
// <= sqrt(lambda_max(M)).
void evaluateMetric
(
    const vectorField& g,
    const tensorField& M,
    const vectorField& w,
    scalarField& G,
    vectorField& a,
    const word& location,
    const vectorField& centres
)
{
    const scalar u = 0.5*std::numeric_limits<scalar>::epsilon();
    const scalar gamma6 = 6*u/(1 - 6*u);

    forAll(g, i)
    {
        scalar Q = g[i] & w[i];

        if (!std::isfinite(Q))
        {
            FatalErrorInFunction
                << "Non-finite grad(activationTime) & "
                << "(M & grad(activationTime)) = " << Q << " at "
                << location << " " << i << ", centre " << centres[i] << ": grad(activationTime) = " << g[i]
                << ", M = " << M[i] << exit(FatalError);
        }

        if (Q < 0)
        {
            const vector magG(cmptMag(g[i]));
            const scalar bound = gamma6*(magG & (cmptMag(M[i]) & magG));

            if (-Q > bound)
            {
                FatalErrorInFunction
                    << "grad(activationTime) & (M & grad(activationTime)) = "
                    << Q << " is negative beyond its rounding error " << bound
                    << " at " << location << " " << i << ", centre "
                    << centres[i] << ": M = " << M[i]
                    << " is not positive definite." << exit(FatalError);
            }

            Q = 0;
        }

        G[i] = Foam::sqrt(Q);
        a[i] = (Q > 0 ? w[i]/G[i] : vector::zero);
    }
}


void evaluateMetric
(
    const volVectorField& g,
    const volTensorField& M,
    const volVectorField& w,
    volScalarField& G,
    volVectorField& a
)
{
    evaluateMetric
    (
        g.primitiveField(),
        M.primitiveField(),
        w.primitiveField(),
        G.primitiveFieldRef(),
        a.primitiveFieldRef(),
        "cell",
        g.mesh().C().primitiveField()
    );

    forAll(G.boundaryField(), patchI)
    {
        evaluateMetric
        (
            g.boundaryField()[patchI],
            M.boundaryField()[patchI],
            w.boundaryField()[patchI],
            G.boundaryFieldRef()[patchI],
            a.boundaryFieldRef()[patchI],
            "face of patch " + G.mesh().boundary()[patchI].name(),
            G.mesh().C().boundaryField()[patchI]
        );
    }
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
    gradActivationTime_
    (
        IOobject
        (
            "grad(" + activationTime_.name() + ")",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedVector("0", dimTime/dimLength, vector::zero)
    ),
    stimulusCellIDs_(0),
    electroProperties_(electroProperties),
    chi_("chi", dimArea/dimVolume, electroProperties),
    Cm_("cm", dimCurrent*dimTime/(dimVoltage*dimArea), electroProperties),
    conductivity_(initialiseConductivity()),
    M_("M", conductivity_/(chi_*Cm_)),
    w_("w", M_ & gradActivationTime_),
    G_("G", sqrt(gradActivationTime_ & w_)),
    c0_("c0", electroProperties),
    a_
    (
        IOobject
        (
            "a",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("a", w_.dimensions()/G_.dimensions(), Zero)
    ),
    u_("u", c0_*a_),
    phiU_("phiU", (fvc::interpolate(u_) & mesh().Sf())),
    divPhiU_("divPhiU", fvc::div(phiU_)),
    eikonalAdvectionDiffusionApproach_
    (
        electroProperties.lookup("eikonalAdvectionDiffusionApproach")
    ),
    useGraphPrePopulation_
    (
        electroProperties.lookupOrDefault<Switch>("useGraphPrePopulation", true)
    ),
    sealedHeartBoundary_
    (
        electroProperties.get<Switch>("sealedHeartBoundary")
    ),
    verificationModelPtr_()
{
    sealedWallPolicy::publish(activationTime_, sealedHeartBoundary_);

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
        checkVerificationBoundaryTypes(activationTime_, verificationModelPtr_());
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

    if (useGraphPrePopulation_) preInitialiseFromSeeds(constrainedCells, constrainedValues);

    const dimensionedScalar one("one", dimless, 1.0);

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
        evaluateMetric(gradActivationTime_, M_, w_, G_, a_);

        if (eikonalAdvectionDiffusionApproach_)
        {
            u_ = c0_*a_;
            phiU_ = (fvc::interpolate(u_) & mesh().Sf());
            divPhiU_ = fvc::div(phiU_);

            fvScalarMatrix activationEqn
            (
               -fvm::laplacian
                (
                    insulatedFaceConductivity
                    (
                        M_, activationTime_, sealedHeartBoundary_
                    )(),
                    activationTime_
                )
              + fvm::div(phiU_, activationTime_)
              + fvm::SuSp(-divPhiU_, activationTime_)
             == one + Smms
            );

            // c0 |grad psi|_M = u & grad psi with u lagged; a bounded div scheme keeps the transport monotone
            activationEqn.relax();
            activationEqn.setValues(constrainedCells, constrainedValues);
            activationEqn.solve("asymmetric_" + activationTime_.name());
        }
        else
        {
            fvScalarMatrix activationEqn
            (
               -fvm::laplacian
                (
                    insulatedFaceConductivity
                    (
                        M_, activationTime_, sealedHeartBoundary_
                    )(),
                    activationTime_
                )
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
