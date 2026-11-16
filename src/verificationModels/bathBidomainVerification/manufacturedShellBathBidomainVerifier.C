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

#include "manufacturedShellBathBidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "fixedValueFvPatchFields.H"
#include "ionicModel.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;
namespace mms = shellBathBidomainMMS;

defineTypeNameAndDebug(manufacturedShellBathBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedShellBathBidomainVerifier,
    dictionary
);

const Enum<manufacturedShellBathBidomainVerifier::outerClosure>
manufacturedShellBathBidomainVerifier::outerClosureNames
({
    { outerClosure::insulated, "insulated" },
    { outerClosure::grounded, "grounded" }
});


namespace
{

const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);

// Tolerance of the geometry checks on the radii of the spheres
const scalar geometryTol = 1e-8;

scalar sampleValue(const mms::SelfTestSample& s, const mms::Params& p)
{
    const point pt(s.x, s.y, s.z);
    switch (s.fn)
    {
        case 0: return mms::Vm(pt, s.t, p);
        case 1: return mms::phiE(pt, s.t, s.zone, p);
        case 2: return mms::phiI(pt, s.t, p);
        case 3: return mms::gradPhiE(pt, s.t, s.zone, p)[s.component];
        case 4: return mms::gradPhiI(pt, s.t, p)[s.component];
        case 5: return mms::gradVm(pt, s.t, p)[s.component];
        case 6: return mms::u1(pt, s.t, p);
        case 7: return mms::u2(pt, s.t, p);
        case 8: return mms::Gi(pt, p)[s.component];
        case 9: return mms::Ge(pt, p)[s.component];
        case 10: return mms::sigmaB(pt, s.zone, p);
        case 11: return mms::divGiGradPhiI(pt, s.t, p);
        case 12: return mms::divGeGradPhiE(pt, s.t, p);
        case 13: return mms::phiESource(pt, s.t, s.zone, p);
        case 14:
        {
            // the flux rows sample the components of Ge grad(phiE)
            vector e(Zero);
            e[s.component] = 1;
            return mms::interfaceFlux(pt, e, s.t, p);
        }
        case 15: return mms::G(pt, p);
    }
    FatalErrorInFunction
        << "Unknown self-test function " << s.fn << exit(FatalError);
    return 0;
}


// Zone of a point by its radius: for the boundary-face values of the exact
// reconstruction, where no cell zone applies (the formulas extend smoothly)
label radialZone(const point& pt, const mms::Params& p)
{
    return mms::zoneOf(pt, p);
}

} // End anonymous namespace


void manufacturedShellBathBidomainVerifier::selfTest()
{
    for (label i = 0; i < mms::nSelfTestSamples; ++i)
    {
        const mms::SelfTestSample& s = mms::selfTestSamples[i];
        mms::Params p =
            s.param == 0 ? mms::defaultParams() : mms::alternateParams();
        p.solution = s.solution;
        const scalar value = sampleValue(s, p);

        if (!(mag(value - s.value) <= 1e-11*max(scalar(1), mag(s.value))))
        {
            FatalErrorInFunction
                << "Bath shell MMS self-test sample " << i
                << " (function " << s.fn << ", parameter set " << s.param
                << ", solution " << s.solution << ", component "
                << s.component << "): " << value << " instead of " << s.value
                << exit(FatalError);
        }
    }

    Info<< "bath shell MMS self-test: " << mms::nSelfTestSamples
        << " samples OK" << endl;
}


manufacturedShellBathBidomainVerifier::manufacturedShellBathBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    p_(),
    closure_(outerClosure::insulated),
    parameterSet_("default"),
    nPerDirection_(0),
    phiEPtr_(nullptr),
    heartCellMapPtr_(nullptr),
    errorsReported_(false),
    beta_(0),
    betaInitialised_(false),
    electrodeNames_(),
    electrodePositions_(),
    zones_()
{
    const dictionary& cfg = dict.subDict("verificationModel");
    for (const entry& e : cfg)
    {
        const word& key = e.keyword();
        if
        (
            key != "type" && key != "N" && key != "outerClosure"
         && key != "parameterSet" && key != "electrodePositions"
        )
        {
            FatalIOErrorInFunction(cfg)
                << type() << " takes the keys N, outerClosure, parameterSet "
                << "and electrodePositions; found " << key << "."
                << exit(FatalIOError);
        }
    }

    const dictionary& electrodes = cfg.subDict("electrodePositions");
    for (const entry& e : electrodes)
    {
        electrodeNames_.append(e.keyword());
        electrodePositions_.append(electrodes.get<point>(e.keyword()));
    }
    if (electrodeNames_.empty())
    {
        FatalIOErrorInFunction(electrodes)
            << "electrodePositions names no electrode." << exit(FatalIOError);
    }

    nPerDirection_ = cfg.get<label>("N");
    if (nPerDirection_ <= 0)
    {
        FatalIOErrorInFunction(cfg)
            << "N must be positive; found " << nPerDirection_ << "."
            << exit(FatalIOError);
    }
    closure_ = outerClosureNames.get("outerClosure", cfg);
    parameterSet_ = cfg.getOrDefault<word>("parameterSet", "default");
    if (parameterSet_ == "default")
    {
        p_ = mms::defaultParams();
    }
    else if (parameterSet_ == "alternate")
    {
        p_ = mms::alternateParams();
    }
    else
    {
        FatalIOErrorInFunction(cfg)
            << "parameterSet must be default or alternate; found "
            << parameterSet_ << "." << exit(FatalIOError);
    }
    p_.solution =
        closure_ == outerClosure::insulated ? mms::insulated : mms::grounded;

    selfTest();
}


void manufacturedShellBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE
)
{
    // The heart-only binding: this solution needs the global phiE of
    // torsoCoupling bath, bound with its heart cell map
    phiEPtr_ = &phiE;
    heartCellMapPtr_ = nullptr;
}


void manufacturedShellBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    phiEPtr_ = &phiE;
    heartCellMapPtr_ = &heartCellMap;
}


void manufacturedShellBathBidomainVerifier::unbindBidomainField()
{
    phiEPtr_ = nullptr;
    heartCellMapPtr_ = nullptr;
}


void manufacturedShellBathBidomainVerifier::initialiseBeta
(
    const ionicModel& model
)
{
    beta_ = requireIonicConstant(typeName, model, "Beta");
    betaInitialised_ = true;
}


const labelList& manufacturedShellBathBidomainVerifier::zones
(
    const fvMesh& baseMesh
) const
{
    if (zones_.size() == baseMesh.nCells())
    {
        return zones_;
    }
    if (!heartCellMapPtr_)
    {
        FatalErrorInFunction
            << typeName << " needs the global phiE of torsoCoupling bath, "
            << "bound with its heart cell map." << exit(FatalError);
    }

    zones_.setSize(baseMesh.nCells());
    zones_ = -1;
    for (const label cellI : *heartCellMapPtr_)
    {
        zones_[cellI] = mms::myocardium;
    }

    const wordList names({"bloodPool", "bath"});
    const label ids[2] = {mms::bloodPool, mms::bath};
    forAll(names, i)
    {
        const label zoneI = baseMesh.cellZones().findZoneID(names[i]);
        if (zoneI < 0)
        {
            FatalErrorInFunction
                << typeName << " needs the cellZone " << names[i] << "."
                << exit(FatalError);
        }
        for (const label cellI : baseMesh.cellZones()[zoneI])
        {
            if (zones_[cellI] >= 0)
            {
                FatalErrorInFunction
                    << "Cell " << cellI << " of cellZone " << names[i]
                    << " is also a heart cell or in another bath zone."
                    << exit(FatalError);
            }
            zones_[cellI] = ids[i];
        }
    }

    // Every cell in a zone, its centre within its own size of the zone's
    // radial range (polyhedral interfaces cut the spheres)
    const vectorField& C = baseMesh.C().primitiveField();
    const scalarField& V = baseMesh.V().field();
    label unassigned = 0;
    label misplaced = 0;
    forAll(zones_, cellI)
    {
        const label z = zones_[cellI];
        if (z < 0)
        {
            ++unassigned;
            continue;
        }
        const scalar r = mag(C[cellI]);
        const scalar rMin =
            z == mms::bloodPool ? 0 : (z == mms::myocardium ? p_.r1 : p_.r2);
        const scalar rMax =
            z == mms::bloodPool ? p_.r1 : (z == mms::myocardium ? p_.r2 : p_.r3);
        const scalar size = Foam::cbrt(V[cellI]);
        if (r < rMin - size || r > rMax + size)
        {
            ++misplaced;
        }
    }
    reduce(unassigned, sumOp<label>());
    reduce(misplaced, sumOp<label>());
    if (unassigned > 0 || misplaced > 0)
    {
        FatalErrorInFunction
            << unassigned << " cells lie in no bath-shell zone and "
            << misplaced << " cells lie farther than their own size outside "
            << "the radial range of their zone." << exit(FatalError);
    }

    return zones_;
}


void manufacturedShellBathBidomainVerifier::validateConfiguration
(
    const ionicModel& model
) const
{
    requireIonicFamily(typeName, model, "monodomainFDAManufactured", 3);

    if (!dict().get<Switch>("sealedHeartBoundary"))
    {
        FatalErrorInFunction
            << typeName << " seals the intracellular space on both spheres "
            << "and needs sealedHeartBoundary true." << exit(FatalError);
    }

}


void manufacturedShellBathBidomainVerifier::validateGeometry
(
    const fvMesh& baseMesh
) const
{
    const labelList& zone = zones(baseMesh);

    // Every point of a cell lies in the radial range of the cell's zone
    const pointField& pts = baseMesh.points();
    const cellList& cells = baseMesh.cells();
    const faceList& faces = baseMesh.faces();
    label nOutside = 0;
    forAll(cells, cellI)
    {
        const label z = zone[cellI];
        const scalar rMin =
            z == mms::bloodPool ? 0 : (z == mms::myocardium ? p_.r1 : p_.r2);
        const scalar rMax =
            z == mms::bloodPool ? p_.r1 : (z == mms::myocardium ? p_.r2 : p_.r3);
        for (const label faceI : cells[cellI])
        {
            for (const label pointI : faces[faceI])
            {
                const scalar r = mag(pts[pointI]);
                if (r < rMin - geometryTol || r > rMax + geometryTol)
                {
                    ++nOutside;
                }
            }
        }
    }
    reduce(nOutside, sumOp<label>());
    if (nOutside > 0)
    {
        FatalErrorInFunction
            << typeName << ": " << nOutside << " cell points lie outside the "
            << "radial range of their zone (bloodPool r <= " << p_.r1
            << ", myocardium " << p_.r1 << " <= r <= " << p_.r2
            << ", bath " << p_.r2 << " <= r <= " << p_.r3 << ")."
            << exit(FatalError);
    }

    // Non-coupled boundary faces: their points on r3
    scalar worstPoint = 0;
    scalar worstCentre = 0;
    for (const polyPatch& pp : baseMesh.boundaryMesh())
    {
        if (pp.coupled() || polyPatch::constraintType(pp.type()))
        {
            continue;
        }
        forAll(pp, i)
        {
            for (const label pointI : pp[i])
            {
                worstPoint = max(worstPoint, mag(mag(pts[pointI]) - p_.r3));
            }
            worstCentre = max(worstCentre, mag(mag(pp.faceCentres()[i]) - p_.r3));
        }
    }
    reduce(worstPoint, maxOp<scalar>());
    reduce(worstCentre, maxOp<scalar>());
    if (worstPoint > geometryTol)
    {
        FatalErrorInFunction
            << typeName << ": a point of the outer boundary lies "
            << worstPoint << " from the sphere r3 = " << p_.r3 << "."
            << exit(FatalError);
    }
    Info<< typeName << ": outer boundary on r3 = " << p_.r3
        << ", largest face-centre distance " << worstCentre << endl;

    // phiE patches against the outer closure
    const volScalarField& phiE = *phiEPtr_;
    forAll(phiE.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = phiE.boundaryField()[patchI];
        if (pf.coupled() || polyPatch::constraintType(pf.patch().patch().type()))
        {
            continue;
        }
        const bool fixedValue = isA<fixedValueFvPatchScalarField>(pf);
        if (closure_ == outerClosure::insulated && fixedValue)
        {
            FatalErrorInFunction
                << typeName << " with outerClosure insulated needs no ground "
                << "patch; phiE is fixedValue on " << pf.patch().name() << "."
                << exit(FatalError);
        }
        if (closure_ == outerClosure::grounded && !fixedValue)
        {
            FatalErrorInFunction
                << typeName << " with outerClosure grounded needs every outer "
                << "patch in groundPatches; " << pf.patch().name()
                << " is " << pf.type() << "." << exit(FatalError);
        }
    }
}


void manufacturedShellBathBidomainVerifier::validateConductivity
(
    const fvMesh& heartMesh,
    const fvMesh& baseMesh
) const
{
    const vectorField& Ch = heartMesh.C().primitiveField();
    const word names[2] =
        {"conductivityIntracellular", "conductivityExtracellular"};
    for (label i = 0; i < 2; ++i)
    {
        const volTensorField& G = heartMesh.lookupObject<volTensorField>(names[i]);
        scalar worst = 0;
        forAll(Ch, cellI)
        {
            const tensor exact
            (
                (i == 0 ? mms::Gi(Ch[cellI], p_) : mms::Ge(Ch[cellI], p_))
              & tensor(I)
            );
            worst = max(worst, mag(G[cellI] - exact)/mag(exact));
        }
        reduce(worst, maxOp<scalar>());
        if (worst > 1e-9)
        {
            FatalErrorInFunction
                << typeName << " needs " << names[i] << " written by "
                << "setManufacturedFields (conductivitySource field); the "
                << "largest relative difference from the manufactured tensor "
                << "is " << worst << "." << exit(FatalError);
        }
    }

    // The bath conductivity of the conductor
    const labelList& zone = zones(baseMesh);
    const volTensorField& sigmaTotal =
        baseMesh.lookupObject<volTensorField>("sigmaTotal");
    const vectorField& C = baseMesh.C().primitiveField();
    scalar worst = 0;
    forAll(C, cellI)
    {
        if (zone[cellI] == mms::myocardium)
        {
            continue;
        }
        const scalar exact = mms::sigmaB(C[cellI], zone[cellI], p_);
        worst = max(worst, mag(sigmaTotal[cellI] - exact*tensor::I)/exact);
    }
    reduce(worst, maxOp<scalar>());
    if (worst > 1e-9)
    {
        FatalErrorInFunction
            << typeName << " needs bodyAndOrgansConductivity written by "
            << "setManufacturedFields; the largest relative difference of the "
            << "bath conductivity from sigma_b is " << worst << "."
            << exit(FatalError);
    }
}


tmp<scalarField> manufacturedShellBathBidomainVerifier::exactPhiE
(
    const fvMesh& baseMesh,
    const scalar t
) const
{
    const labelList& zone = zones(baseMesh);
    const vectorField& C = baseMesh.C().primitiveField();
    tmp<scalarField> tphi(new scalarField(C.size()));
    scalarField& phi = tphi.ref();
    forAll(C, cellI)
    {
        phi[cellI] = mms::phiE(C[cellI], t, zone[cellI], p_);
    }
    return tphi;
}


wordList manufacturedShellBathBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedShellBathBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedShellBathBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedShellBathBidomainVerifier::preProcess
(
    ionicModel& model,
    volScalarField& Vm,
    PtrList<volScalarField>& fields
)
{
    const wordList names = preProcessFieldNames(model);
    validateConfiguration(model);
    initialiseBeta(model);

    // The myocardium calls preProcess with its own phiE first and again once
    // the bath binds the global phiE with the heart cell map; the bath checks
    // run on the second call
    const fvMesh& heartMesh = Vm.mesh();
    if (heartCellMapPtr_)
    {
        const fvMesh& baseMesh = phiEPtr_->mesh();
        validateGeometry(baseMesh);
        validateConductivity(heartMesh, baseMesh);
    }

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    const vectorField& C = heartMesh.C().primitiveField();
    const scalar t = heartMesh.time().value();

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = fields[u1Field].primitiveFieldRef();
    scalarField& u2I = fields[u2Field].primitiveFieldRef();
    scalarField& u3I = fields[u3Field].primitiveFieldRef();
    forAll(C, cellI)
    {
        VmI[cellI] = mms::Vm(C[cellI], t, p_);
        u1I[cellI] = mms::u1(C[cellI], t, p_);
        u2I[cellI] = mms::u2(C[cellI], t, p_);
        u3I[cellI] = 0;
    }
    Vm.correctBoundaryConditions();
    fields[u1Field].correctBoundaryConditions();
    fields[u2Field].correctBoundaryConditions();
    fields[u3Field].correctBoundaryConditions();

    if (phiEPtr_ && heartCellMapPtr_)
    {
        phiEPtr_->primitiveFieldRef() = exactPhiE(phiEPtr_->mesh(), t);
        phiEPtr_->correctBoundaryConditions();
    }
    else if (phiEPtr_)
    {
        const vectorField& Cp = phiEPtr_->mesh().C().primitiveField();
        scalarField& phiEI = phiEPtr_->primitiveFieldRef();
        forAll(Cp, cellI)
        {
            phiEI[cellI] = mms::phiE(Cp[cellI], t, mms::myocardium, p_);
        }
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedShellBathBidomainVerifier::addManufacturedPdeSource
(
    volScalarField& sourceField,
    const volTensorField& conductivity,
    const scalar evaluationTime
) const
{
    if (!betaInitialised_)
    {
        FatalErrorInFunction
            << "The manufactured source was requested before preProcess "
            << "read the ionic constant Beta." << exit(FatalError);
    }

    const vectorField& C = sourceField.mesh().C().primitiveField();
    scalarField& source = sourceField.primitiveFieldRef();

    // The FDA ionic model contributes -Cm dV/dt + beta V/chi, so the source
    // is beta Vm - div(Gi grad phiI)
    forAll(C, cellI)
    {
        source[cellI] +=
            beta_*mms::Vm(C[cellI], evaluationTime, p_)
          - mms::divGiGradPhiI(C[cellI], evaluationTime, p_);
    }

    sourceField.correctBoundaryConditions();
}


void manufacturedShellBathBidomainVerifier::addManufacturedPhiESource
(
    volScalarField& rhs,
    const scalar evaluationTime
) const
{
    const fvMesh& baseMesh = rhs.mesh();
    const labelList& zone = zones(baseMesh);
    const vectorField& C = baseMesh.C().primitiveField();
    const scalarField& V = baseMesh.V().field();

    scalarField source(C.size());
    forAll(C, cellI)
    {
        source[cellI] =
            mms::phiESource(C[cellI], evaluationTime, zone[cellI], p_);
    }

    // Insulated: the continuum source integrates to zero over heart and bath;
    // its cell-centre samples do not, by a quadrature error. Removing their
    // volume-weighted mean makes the pure-Neumann system compatible on the
    // mesh. Grounded: no compatibility condition applies, and a shift of the
    // source would change the solution
    if (closure_ == outerClosure::insulated)
    {
        source -= gSum(source*V)/gSum(V);
    }

    rhs.primitiveFieldRef() += source;
}


void manufacturedShellBathBidomainVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " writes the manufactured conductivities of the "
            << "undecomposed case; run it in serial." << exit(FatalError);
    }

    const unsigned int oldPrecision = IOstream::defaultPrecision(17);
    const vectorField& C = mesh.C().primitiveField();

    const auto io = [&mesh](const word& name)
    {
        return IOobject(name, "0", mesh, IOobject::NO_READ, IOobject::NO_WRITE);
    };
    const auto copyInternalToPatches = [&mesh](auto& field)
    {
        forAll(mesh.boundary(), patchI)
        {
            if (polyPatch::constraintType(mesh.boundaryMesh()[patchI].type()))
            {
                continue;
            }
            field.boundaryFieldRef()[patchI] ==
                field.boundaryField()[patchI].patchInternalField();
        }
    };

    // The heart tensors everywhere: a heart cell whose centre lies beyond a
    // polyhedral interface carries the tensor its sources use
    for (label i = 0; i < 2; ++i)
    {
        volSymmTensorField G
        (
            io(i == 0 ? "ConductivityIntracellular" : "ConductivityExtracellular"),
            mesh,
            dimensionedSymmTensor(conductivityDim, Zero),
            "calculated"
        );
        forAll(C, cellI)
        {
            G[cellI] =
                i == 0 ? mms::Gi(C[cellI], p_) : mms::Ge(C[cellI], p_);
        }
        copyInternalToPatches(G);
        G.write();
        Info<< typeName << ": wrote " << G.name() << " in 0." << endl;
    }

    // The bath conductivity by cellZone; the heart cells take the outer
    // bath's formula, which the conductor does not use
    labelList zone(mesh.nCells(), mms::bath);
    const label poolI = mesh.cellZones().findZoneID("bloodPool");
    if (poolI < 0)
    {
        FatalErrorInFunction
            << typeName << " needs the cellZone bloodPool." << exit(FatalError);
    }
    for (const label cellI : mesh.cellZones()[poolI])
    {
        zone[cellI] = mms::bloodPool;
    }
    volScalarField sigma
    (
        io("bodyAndOrgansConductivity"),
        mesh,
        dimensionedScalar(conductivityDim, Zero),
        "calculated"
    );
    forAll(C, cellI)
    {
        sigma[cellI] = mms::sigmaB(C[cellI], zone[cellI], p_);
    }
    copyInternalToPatches(sigma);
    sigma.write();
    Info<< typeName << ": wrote " << sigma.name() << " in 0." << endl;

    IOstream::defaultPrecision(oldPrecision);
}


void manufacturedShellBathBidomainVerifier::postProcess
(
    const ionicModel& model,
    const volScalarField& Vm,
    const PtrList<volScalarField>& fields
)
{
    if (!shouldPostProcess(model, Vm))
    {
        return;
    }
    validateConfiguration(model);
    if (!phiEPtr_ || !heartCellMapPtr_)
    {
        FatalErrorInFunction
            << typeName << " needs the global phiE of torsoCoupling bath, "
            << "bound with its heart cell map." << exit(FatalError);
    }

    const wordList names = requiredPostProcessFieldNames(model);
    const label u1Idx = requireFieldIndex(names, "u1", "postProcess");
    const label u2Idx = requireFieldIndex(names, "u2", "postProcess");

    const fvMesh& heartMesh = Vm.mesh();
    const fvMesh& baseMesh = phiEPtr_->mesh();
    const labelUList& heartCellMap = *heartCellMapPtr_;
    const labelList& zone = zones(baseMesh);
    const Time& time = heartMesh.time();
    const scalar t = time.value();
    const vectorField& Ch = heartMesh.C().primitiveField();
    const vectorField& Cb = baseMesh.C().primitiveField();

    // Heart fields
    scalarField VmExact(Ch.size()), u1Exact(Ch.size()), u2Exact(Ch.size());
    scalarField phiIExact(Ch.size());
    vectorField gradVmExact(Ch.size()), gradPhiIExact(Ch.size());
    forAll(Ch, cellI)
    {
        VmExact[cellI] = mms::Vm(Ch[cellI], t, p_);
        u1Exact[cellI] = mms::u1(Ch[cellI], t, p_);
        u2Exact[cellI] = mms::u2(Ch[cellI], t, p_);
        phiIExact[cellI] = mms::phiI(Ch[cellI], t, p_);
        gradVmExact[cellI] = mms::gradVm(Ch[cellI], t, p_);
        gradPhiIExact[cellI] = mms::gradPhiI(Ch[cellI], t, p_);
    }

    // Global phiE: compared as it is when grounded; after removing the
    // volume mean of its error when insulated (the pinned reference cell
    // fixes an arbitrary constant)
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    scalarField phiEExact(exactPhiE(baseMesh, t));
    const scalar gaugeOffset =
        closure_ == outerClosure::insulated
      ? computeVolumeMeanError(baseMesh, phiEValues, phiEExact)
      : 0;
    phiEExact += gaugeOffset;
    const word phiERow
    (
        closure_ == outerClosure::insulated ? "phiE_gauge" : "phiE"
    );

    vectorField gradPhiEExact(Cb.size());
    boolList isHeart(Cb.size(), false);
    forAll(Cb, cellI)
    {
        gradPhiEExact[cellI] = mms::gradPhiE(Cb[cellI], t, zone[cellI], p_);
        isHeart[cellI] = zone[cellI] == mms::myocardium;
    }

    // phiI = Vm + phiE on the heart cells, with the same gauge
    const scalarField& VmValues = Vm.primitiveField();
    scalarField phiIValues(VmValues.size());
    forAll(phiIValues, i)
    {
        phiIValues[i] = VmValues[i] + phiEValues[heartCellMap[i]];
    }
    phiIExact += gaugeOffset;

    const auto VmNorms = computeNorms(heartMesh, VmValues, VmExact);
    const auto u1Norms =
        computeNorms(heartMesh, fields[u1Idx].primitiveField(), u1Exact);
    const auto u2Norms =
        computeNorms(heartMesh, fields[u2Idx].primitiveField(), u2Exact);
    const auto phiINorms = computeNorms(heartMesh, phiIValues, phiIExact);
    const auto phiENorms = computeNorms(baseMesh, phiEValues, phiEExact);

    // phiE by zone, each normalised by the total volume
    List<Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>> rows;
    const word zoneNames[3] = {"pool", "heart", "bath"};
    for (label z = 0; z < 3; ++z)
    {
        boolList mask(Cb.size());
        forAll(mask, i)
        {
            mask[i] = zone[i] == z;
        }
        rows.append
        (
            Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>
            (
                phiERow + "_" + zoneNames[z],
                computeMaskedNorms(baseMesh, phiEValues, phiEExact, mask)
            )
        );
    }

    // Bulk, boundary and interface-band rows
    const tmp<volVectorField> tgradVm(fvc::grad(Vm));
    const tmp<volVectorField> tgradReconVm
    (
        exactFieldGradient
        (
            heartMesh,
            Vm.name(),
            [&](const point& c) { return mms::Vm(c, t, p_); }
        )
    );
    rows.append
    (
        bulkBoundaryRows
        (
            heartMesh, "Vm", VmValues, VmExact, tgradVm().primitiveField(),
            tgradReconVm().primitiveField(), gradVmExact
        )
    );

    const volScalarField& heartPhiE =
        heartMesh.lookupObject<volScalarField>("phiE");
    const vectorField gradPhiI
    (
        tgradVm().primitiveField() + fvc::grad(heartPhiE)().primitiveField()
    );
    const tmp<volVectorField> tgradReconPhiI
    (
        exactFieldGradient
        (
            heartMesh,
            "phiI",
            [&](const point& c) { return mms::phiI(c, t, p_) + gaugeOffset; }
        )
    );
    rows.append
    (
        bulkBoundaryRows
        (
            heartMesh, "phiI", phiIValues, phiIExact, gradPhiI,
            tgradReconPhiI().primitiveField(), gradPhiIExact
        )
    );

    const tmp<volVectorField> tgradPhiE(consumedPhiEGradient(*phiEPtr_));
    const tmp<volVectorField> tgradReconPhiE
    (
        exactFieldGradient
        (
            baseMesh,
            phiEPtr_->name(),
            [&](const point& c)
            {
                return mms::phiE(c, t, radialZone(c, p_), p_) + gaugeOffset;
            }
        )
    );
    rows.append
    (
        globalPotentialRows
        (
            baseMesh, isHeart, phiERow, phiEValues, phiEExact,
            tgradPhiE().primitiveField(), tgradReconPhiE().primitiveField(),
            gradPhiEExact
        )
    );

    // Electrodes: the solved phiE of the cell that holds each point against
    // the exact phiE at the cell centre and at the point, gathered from the
    // rank that owns the cell
    const label nE = electrodeNames_.size();
    scalarField eComputed(nE, Zero), eCentroid(nE, Zero), ePoint(nE, Zero);
    labelList eFound(nE, Zero);
    forAll(electrodePositions_, e)
    {
        const point& pt = electrodePositions_[e];
        const label cellI = baseMesh.findCell(pt);
        if (cellI < 0)
        {
            continue;
        }
        eFound[e] = 1;
        eComputed[e] = phiEValues[cellI];
        eCentroid[e] = phiEExact[cellI];
        ePoint[e] = mms::phiE(pt, t, zone[cellI], p_) + gaugeOffset;
    }
    Pstream::listCombineReduce(eFound, plusEqOp<label>());
    Pstream::listCombineReduce(eComputed, plusEqOp<scalar>());
    Pstream::listCombineReduce(eCentroid, plusEqOp<scalar>());
    Pstream::listCombineReduce(ePoint, plusEqOp<scalar>());
    forAll(eFound, e)
    {
        if (eFound[e] != 1)
        {
            FatalErrorInFunction
                << typeName << ": electrode " << electrodeNames_[e] << " at "
                << electrodePositions_[e] << " lies in " << eFound[e]
                << " cells of the mesh." << exit(FatalError);
        }
    }

    // The body surface: the electrodes named surface*, as one row (mean,
    // RMS and largest error against the centroid value), and one row per
    // electrode
    scalar sAbs = 0, sSq = 0, sMax = 0;
    label nSurface = 0;
    forAll(electrodeNames_, e)
    {
        const scalar err = mag(eComputed[e] - eCentroid[e]);
        rows.append
        (
            Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>
            (
                phiERow + "_electrode_" + electrodeNames_[e],
                Tuple2<Tuple2<scalar, scalar>, scalar>
                (
                    Tuple2<scalar, scalar>(err, err), err
                )
            )
        );
        if (electrodeNames_[e].starts_with("surface"))
        {
            sAbs += err;
            sSq += err*err;
            sMax = max(sMax, err);
            ++nSurface;
        }
    }
    if (nSurface > 0)
    {
        rows.append
        (
            Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>
            (
                phiERow + "_bodySurface",
                Tuple2<Tuple2<scalar, scalar>, scalar>
                (
                    Tuple2<scalar, scalar>
                    (
                        sAbs/nSurface, Foam::sqrt(sSq/nSurface)
                    ),
                    sMax
                )
            )
        );
    }

    const label nCells = returnReduce(baseMesh.nCells(), sumOp<label>());
    const scalar h = Foam::cbrt(gSum(baseMesh.V())/max(nCells, label(1)));
    const word cellKind(meshCellKind(baseMesh));
    const fileName outputDir(postProcessingDir(time));
    const fileName outputFile
    (
        outputDir/studyFileName("3D_" + cellKind, nPerDirection_)
    );

    if (Pstream::master())
    {
        mkDir(outputDir);
        OFstream os(outputFile);
        os.precision(12);
        os  << "# Bath-shell bidomain manufactured solution error summary\n"
            << "# mesh " << cellKind << "\n"
            << "# cellsPerDirection " << nPerDirection_ << "\n"
            << "# h " << h << "\n"
            << "# cells " << nCells << "\n"
            << "# steps " << max(label(0), time.timeIndex()) << "\n"
            << "# dt " << time.deltaTValue() << "\n"
            << "# time " << t << "\n"
            << "# Beta " << beta_ << "\n"
            << "# outerClosure " << outerClosureNames[closure_] << "\n"
            << "# parameterSet " << parameterSet_ << "\n"
            << "# field L1 L2 Linf\n";
        os  << "# electrode x y z computed exactCentroid exactPoint\n";
        forAll(electrodeNames_, e)
        {
            const point& pt = electrodePositions_[e];
            os  << "# electrode " << electrodeNames_[e] << " " << pt.x() << " "
                << pt.y() << " " << pt.z() << " " << eComputed[e] << " "
                << eCentroid[e] << " " << ePoint[e] << "\n";
        }
        writeNormRow(os, "Vm", VmNorms);
        writeNormRow(os, "u1", u1Norms);
        writeNormRow(os, "u2", u2Norms);
        writeNormRow(os, "phiI", phiINorms);
        writeNormRow(os, phiERow, phiENorms);
        for (const auto& row : rows)
        {
            writeNormRow(os, row.first(), row.second());
        }

        errorTable table;
        table.titlePrefix =
            "Bath-shell bidomain manufactured-solution error summary (t = ";
        table.t = t;
        table.addRow("Vm        ", VmNorms);
        table.addRow("u1        ", u1Norms);
        table.addRow("u2        ", u2Norms);
        table.addRow("phiI      ", phiINorms);
        table.addRow(closure_ == outerClosure::insulated ? "phiE_gauge" : "phiE      ", phiENorms);
        table.splitRows = rows;
        writeErrorTableConsole(Info, table);
        Info<< "Written to " << outputFile << endl;
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
