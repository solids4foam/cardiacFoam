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

#include "manufacturedLeadFieldECGVerifier.H"

#include "addToRunTimeSelectionTable.H"
#include "ecgModelIO.H"
#include "fvMatrices.H"
#include "emptyFvPatch.H"
#include "PstreamReduceOps.H"
#include <limits>

namespace Foam
{

defineTypeNameAndDebug(manufacturedLeadFieldECGVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    manufacturedLeadFieldECGVerifier,
    dictionary
);


const Enum<manufacturedLeadFieldECGVerifier::referenceType>
manufacturedLeadFieldECGVerifier::referenceTypeNames
({
    { referenceType::slab, "slab" },
    { referenceType::layered, "layered" },
    { referenceType::sphericalShell, "sphericalShell" }
});


const wordList manufacturedLeadFieldECGVerifier::quantityNames_
({
    "leadField-direct",
    "leadFieldStatic-directStatic",
    "leadFieldStatic-exactPoint",
    "leadFieldStatic-exactCentroid",
    "directStatic-exactPoint",
    "directStatic-exactCentroid",
    "localisation"
});


const wordList manufacturedLeadFieldECGVerifier::fieldRegionNames_
({
    "phiStatic",
    "phiStatic_heart",
    "phiStatic_torso",
    "phiStatic_interfaceBand",
    "phiStatic_surface"
});


namespace
{


// Layered: tolerance of the face normals and the cell extents (unit box)
const scalar layeredGeometryTol = 1e-9;

// sphericalShell: relative tolerance of the case conductivities against the
// solution's at the cell centres
const scalar shellConductivityTol = 1e-10;

} // End anonymous namespace


manufacturedLeadFieldECGVerifier::manufacturedLeadFieldECGVerifier
(
    const electroStateProvider& stateProvider,
    const word& domainName,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    ecgVerificationModel
    (
        stateProvider, domainName, electrodeNames, electrodePositions
    ),
    reference_(referenceType::slab),
    leadFieldPtr_(),
    slabK_(0),
    layers_(),
    alpha_(0),
    shellParams_(shellTransmissionMMS::defaultParams()),
    leadFields_(),
    shellZones_(),
    electrodeCells_(),
    electrodeCentroids_(),
    referenceCentroid_(point::zero),
    directPhi_(),
    directStaticPhi_(),
    heartPhi_(),
    outputPtr_(),
    sampleCount_(0),
    exactScale_(0),
    electrodeNorms_
    (
        quantityNames_.size(),
        List<errorNorms>(electrodePositions.size())
    ),
    maxNorms_(quantityNames_.size()),
    relativeGapNorms_(),
    fieldErrors_(fieldRegionNames_.size(), FixedList<scalar, 3>(Zero)),
    summaryWritten_(false)
{
    read(dict);
}


ecgVerificationModel::Requirements
manufacturedLeadFieldECGVerifier::requirements() const
{
    Requirements needs;
    needs.needVm = true;
    return needs;
}


bool manufacturedLeadFieldECGVerifier::read(const dictionary& dict)
{
    const dictionary& verifierDict = dict.subDict("verificationModel");
    const referenceType selected =
        referenceTypeNames.get("reference", verifierDict);

    if (leadFieldPtr_.valid() && selected != reference_)
    {
        FatalErrorInFunction
            << "The reference of " << type() << " is fixed at the first "
            << "sample." << exit(FatalError);
    }
    reference_ = selected;

    if (reference_ == referenceType::layered)
    {
        const FixedList<scalar, 2> x
        (
            verifierDict.get<FixedList<scalar, 2>>("layerInterfaces")
        );
        if (!(x[0] < 0) || !(x[1] > 1))
        {
            FatalIOErrorInFunction(verifierDict)
                << "layerInterfaces " << x << " must be (xL xR) with xL < 0 "
                << "and xR > 1." << exit(FatalIOError);
        }
        layers_.xLeft = x[0];
        layers_.xRight = x[1];
    }
    else if (verifierDict.found("layerInterfaces"))
    {
        FatalIOErrorInFunction(verifierDict)
            << "layerInterfaces applies to reference layered only."
            << exit(FatalIOError);
    }

    return true;
}


void manufacturedLeadFieldECGVerifier::initialiseLayered
(
    const leadFieldConductor& lfc
)
{
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    if (!lfc.groundPatches().empty())
    {
        FatalErrorInFunction
            << "The layered reference is pinned: it needs no ground patches, "
            << "but has " << lfc.groundPatches() << "." << exit(FatalError);
    }

    // alpha = I n_x on every surface-current face; no net current
    const wordList& patches = lfc.surfaceCurrentPatches();
    const scalarList& currents = lfc.surfaceCurrentValues();
    scalar alphaMin = GREAT;
    scalar alphaMax = -GREAT;
    scalar netCurrent = 0;
    scalar grossCurrent = 0;
    label nFaces = 0;
    forAll(patches, i)
    {
        const fvPatch& patch = baseMesh.boundary()[patches[i]];
        const vectorField nf(patch.nf());
        forAll(nf, faceI)
        {
            if (mag(mag(nf[faceI].x()) - 1) > layeredGeometryTol)
            {
                FatalErrorInFunction
                    << "The layered reference needs surface-current faces "
                    << "normal to x; patch '" << patches[i] << "' has the "
                    << "face normal " << nf[faceI] << "." << exit(FatalError);
            }
            const scalar a = currents[i]*nf[faceI].x();
            alphaMin = min(alphaMin, a);
            alphaMax = max(alphaMax, a);
            netCurrent += currents[i]*patch.magSf()[faceI];
            grossCurrent += mag(currents[i])*patch.magSf()[faceI];
            ++nFaces;
        }
    }
    reduce(alphaMin, minOp<scalar>());
    reduce(alphaMax, maxOp<scalar>());
    reduce(netCurrent, sumOp<scalar>());
    reduce(grossCurrent, sumOp<scalar>());
    reduce(nFaces, sumOp<label>());

    if (nFaces == 0)
    {
        alpha_ = 0;
    }
    else
    {
        if (alphaMax - alphaMin > 1e-12*max(mag(alphaMin), mag(alphaMax)))
        {
            FatalErrorInFunction
                << "The layered reference needs current times n_x equal on "
                << "every surface-current face; it ranges from " << alphaMin
                << " to " << alphaMax << " on " << patches << "."
                << exit(FatalError);
        }
        if (mag(netCurrent) > 1e-12*grossCurrent)
        {
            FatalErrorInFunction
                << "The layered reference is insulated apart from the "
                << "current along x: the net surface current " << netCurrent
                << " on " << patches << " must be zero." << exit(FatalError);
        }
        alpha_ = 0.5*(alphaMin + alphaMax);
    }

    // With a current the exact heart slope alpha/(Gi + Ge)_xx at the walls
    // carries the intracellular flux alpha Gi_xx/(Gi + Ge)_xx, which a
    // sealed Vm admits and an insulated bidomain intracellular potential
    // does not
    if (alpha_ != 0 && !lfc.conductor().totalWallConductivity())
    {
        FatalErrorInFunction
            << "The layered reference with a surface current (alpha "
            << alpha_ << ") is a monodomain solution; the bidomain conductor "
            << "insulates the intracellular potential at the heart wall."
            << exit(FatalError);
    }

    // (Gi + Ge)_xx; slabK_ has checked Gi and Ge
    layers_.sigmaHeart =
        verificationUtils::uniformTensor
        (
            lfc.Gi().primitiveField() + lfc.Ge().primitiveField(),
            "layered reference"
        ).xx();

    // Each torso cell lies in one layer; one isotropic value per layer
    const passiveConductor& conductor = lfc.conductor();
    boolList heart(baseMesh.nCells(), false);
    for (const label cellI : conductor.heartCellToBaseCell())
    {
        heart[cellI] = true;
    }

    const scalar xL = layers_.xLeft;
    const scalar xR = layers_.xRight;
    const scalar tol = layeredGeometryTol;
    const tensorField& sigma = conductor.sigmaTotal().primitiveField();
    const pointField& points = baseMesh.points();
    List<DynamicList<tensor>> layerSigma(4);

    forAll(sigma, cellI)
    {
        if (heart[cellI])
        {
            continue;
        }

        scalar xMin = GREAT;
        scalar xMax = -GREAT;
        for (const label pointI : baseMesh.cellPoints()[cellI])
        {
            xMin = min(xMin, points[pointI].x());
            xMax = max(xMax, points[pointI].x());
        }

        label layer = -1;
        if (xMax <= xL + tol)
        {
            layer = 0;
        }
        else if (xMin >= xL - tol && xMax <= tol)
        {
            layer = 1;
        }
        else if (xMin >= 1 - tol && xMax <= xR + tol)
        {
            layer = 2;
        }
        else if (xMin >= xR - tol)
        {
            layer = 3;
        }
        else
        {
            FatalErrorInFunction
                << "Torso cell " << cellI << " spans " << xMin << " <= x <= "
                << xMax << ", across a layer interface of the layered "
                << "reference (" << xL << ", 0, 1, " << xR << ")."
                << exit(FatalError);
        }
        layerSigma[layer].append(sigma[cellI]);
    }

    const wordList layerNames
    ({
        "outer left", "inner left", "inner right", "outer right"
    });
    scalarList layerValue(4);
    forAll(layerSigma, layer)
    {
        const tensor s = verificationUtils::uniformTensor
        (
            tensorField(layerSigma[layer]),
            "layered reference, " + layerNames[layer] + " layer"
        );
        if (mag(s - s.xx()*tensor::I) > 1e-10*mag(s))
        {
            FatalErrorInFunction
                << "The layered reference needs an isotropic torso "
                << "conductivity; the " << layerNames[layer] << " layer has "
                << s << "." << exit(FatalError);
        }
        layerValue[layer] = s.xx();
    }
    layers_.sigmaOuterLeft = layerValue[0];
    layers_.sigmaInnerLeft = layerValue[1];
    layers_.sigmaInnerRight = layerValue[2];
    layers_.sigmaOuterRight = layerValue[3];
}


void manufacturedLeadFieldECGVerifier::initialiseShell
(
    const leadFieldConductor& lfc
)
{
    using namespace shellTransmissionMMS;

    const passiveConductor& conductor = lfc.conductor();
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    if
    (
        !lfc.groundPatches().empty()
     || !lfc.surfaceCurrentPatches().empty()
    )
    {
        FatalErrorInFunction
            << "The sphericalShell reference is pinned and takes its outer "
            << "current from the solution: it needs no ground or "
            << "surface-current patches, but has " << lfc.groundPatches()
            << " and " << lfc.surfaceCurrentPatches() << "."
            << exit(FatalError);
    }
    if (conductor.totalWallConductivity())
    {
        FatalErrorInFunction
            << "The sphericalShell reference seals the intracellular "
            << "potential at the heart wall: it needs the bidomain conductor."
            << exit(FatalError);
    }

    // Zones from the cellZones against the conductor's heart cells, and the
    // torso conductivity against sigmaB
    shellZones_ = manufacturedLeadFieldShellZones(baseMesh);
    boolList heart(baseMesh.nCells(), false);
    for (const label cellI : conductor.heartCellToBaseCell())
    {
        heart[cellI] = true;
    }
    const vectorField& C = baseMesh.C().primitiveField();
    const tensorField& sigma = conductor.sigmaTotal().primitiveField();
    label misplaced = 0;
    scalar worst = 0;
    forAll(C, cellI)
    {
        const label zone = shellZones_[cellI];
        if (heart[cellI] != (zone == myocardium))
        {
            ++misplaced;
        }
        else if (!heart[cellI])
        {
            const tensor exact = sigmaB(C[cellI], zone, shellParams_)*tensor::I;
            worst = max(worst, mag(sigma[cellI] - exact)/mag(exact));
        }
    }

    // Gi and Ge of the heart cells
    const vectorField& Ch = lfc.Gi().mesh().C().primitiveField();
    const tensorField& GiI = lfc.Gi().primitiveField();
    const tensorField& GeI = lfc.Ge().primitiveField();
    forAll(Ch, cellI)
    {
        const tensor gi = Gi(Ch[cellI], shellParams_) & tensor::I;
        const tensor ge = Ge(Ch[cellI], shellParams_) & tensor::I;
        worst = max(worst, mag(GiI[cellI] - gi)/mag(gi));
        worst = max(worst, mag(GeI[cellI] - ge)/mag(ge));
    }
    reduce(misplaced, sumOp<label>());
    reduce(worst, maxOp<scalar>());

    if (misplaced > 0)
    {
        FatalErrorInFunction
            << misplaced << " cells of the cellZone myocardium differ from "
            << "the conductor's heart cells." << exit(FatalError);
    }
    if (worst > shellConductivityTol)
    {
        FatalErrorInFunction
            << "The case conductivities differ from the sphericalShell "
            << "solution's at the cell centres by up to " << worst
            << " (relative); write them with setManufacturedFields."
            << exit(FatalError);
    }

    forAll(electrodePositions_, e)
    {
        if (zoneOf(electrodePositions_[e], shellParams_) == myocardium)
        {
            FatalErrorInFunction
                << "Electrode " << electrodeNames_[e] << " at "
                << electrodePositions_[e] << " lies in the myocardium of the "
                << "sphericalShell solution." << exit(FatalError);
        }
    }

    leadFields_.setSize(electrodeNames_.size());
    forAll(electrodeNames_, e)
    {
        leadFields_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "leadField_" + electrodeNames_[e],
                    baseMesh.time().constant(),
                    baseMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                baseMesh
            )
        );
    }
}


tmp<volScalarField> manufacturedLeadFieldECGVerifier::shellSource
(
    const scalar t,
    const dimensionSet& dims
) const
{
    using namespace shellTransmissionMMS;

    const fvMesh& baseMesh = stateProvider_.baseMesh();

    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "shellPhiESource",
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            baseMesh,
            dimensionedScalar(dims, 0),
            "zeroGradient"
        )
    );
    scalarField& S = tS.ref().primitiveFieldRef();

    // div(K grad phi) = b + S in every zone
    const vectorField& C = baseMesh.C().primitiveField();
    forAll(C, cellI)
    {
        S[cellI] = phiESource(C[cellI], t, shellZones_[cellI], shellParams_);
    }

    // The exact outward current g on the insulated outer boundary moves into
    // the owner cell's balance: S -= g |Sf|/V
    const scalarField& V = baseMesh.V().field();
    forAll(baseMesh.boundary(), patchI)
    {
        const fvPatch& patch = baseMesh.boundary()[patchI];
        if (patch.coupled() || isA<emptyFvPatch>(patch))
        {
            continue;
        }
        const vectorField nf(patch.nf());
        const vectorField& Cf = patch.Cf();
        const scalarField& magSf = patch.magSf();
        const labelUList& cells = patch.faceCells();
        forAll(cells, faceI)
        {
            const point& pf = Cf[faceI];
            const scalar g =
                sigmaB(pf, bath, shellParams_)
               *(nf[faceI] & gradPhiE(pf, t, bath, shellParams_));
            S[cells[faceI]] -= g*magSf[faceI]/V[cells[faceI]];
        }
    }

    // The pure-Neumann problem needs a source of zero integral
    S -= gSum(V*S)/gSum(V);

    return tS;
}


void manufacturedLeadFieldECGVerifier::initialise()
{
    // The myocardium solver Coeffs of this case
    const IOdictionary& electroProperties =
        mesh_.time().lookupObject<IOdictionary>("electroProperties");
    const dictionary& electroCoeffs = electroProperties.subDict
    (
        electroProperties.get<word>("myocardiumSolver") + "Coeffs"
    );

    leadFieldPtr_.reset(new leadFieldECG(electroCoeffs, domainName_));
    leadFieldPtr_->load(stateProvider_, electrodeNames_);

    const leadFieldConductor& lfc = leadFieldPtr_->conductor();
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    switch (reference_)
    {
        case referenceType::slab:
        {
            if
            (
                !lfc.groundPatches().empty()
             || !lfc.surfaceCurrentPatches().empty()
            )
            {
                FatalErrorInFunction
                    << "The slab reference is insulated: it needs no ground "
                    << "or surface-current patches, but has "
                    << lfc.groundPatches() << " and "
                    << lfc.surfaceCurrentPatches() << "."
                    << exit(FatalError);
            }
            slabK_ = manufacturedLeadFieldSlabK(lfc.Gi(), lfc.Ge());
            break;
        }
        case referenceType::layered:
        {
            slabK_ = manufacturedLeadFieldSlabK(lfc.Gi(), lfc.Ge());
            initialiseLayered(lfc);
            break;
        }
        case referenceType::sphericalShell:
        {
            initialiseShell(lfc);
            break;
        }
    }

    // Electrode cells, their centroids and the centroid of the pinned cell
    electrodeCells_ = lfc.electrodeCells();
    electrodeCentroids_.setSize(electrodeCells_.size(), point::zero);
    forAll(electrodeCells_, e)
    {
        if (electrodeCells_[e] >= 0)
        {
            electrodeCentroids_[e] = baseMesh.C()[electrodeCells_[e]];
        }
    }
    for (point& c : electrodeCentroids_)
    {
        reduce(c, sumOp<point>());
    }

    if (lfc.groundPatches().empty())
    {
        const label refCell = lfc.referenceCell();
        referenceCentroid_ =
            refCell >= 0 ? baseMesh.C()[refCell] : point::zero;
        reduce(referenceCentroid_, sumOp<point>());
    }

    directPhi_.reset(lfc.forwardPotential().ptr());
    directStaticPhi_.reset(lfc.forwardPotential().ptr());
    heartPhi_.reset(lfc.heartPotential().ptr());

    if (Pstream::master())
    {
        wordList columns;
        for (const word& name : electrodeNames_)
        {
            for
            (
                const word& column
              : {
                    "leadField", "direct", "exactPoint", "exactCentroid",
                    "leadFieldStatic", "directStatic"
                }
            )
            {
                columns.append(word(column) + "_" + name);
            }
        }
        outputPtr_ = ecgModelIO::openTimeSeries
        (
            mesh_.time().globalPath()/"postProcessing",
            "manufacturedLeadFieldECG_" + domainName_ + ".dat",
            columns
        );
        outputPtr_->precision(std::numeric_limits<scalar>::max_digits10);
    }

    Info<< type() << ": reference " << referenceTypeNames[reference_];
    if (reference_ != referenceType::sphericalShell)
    {
        Info<< ", k " << slabK_;
    }
    if (reference_ == referenceType::layered)
    {
        Info<< ", alpha " << alpha_ << ", interfaces " << layers_.xLeft
            << " " << layers_.xRight << ", sigma outer left "
            << layers_.sigmaOuterLeft << ", inner left "
            << layers_.sigmaInnerLeft << ", heart " << layers_.sigmaHeart
            << ", inner right " << layers_.sigmaInnerRight
            << ", outer right " << layers_.sigmaOuterRight;
    }
    Info<< ", electrode cell centroids " << electrodeCentroids_ << endl;
}


scalar manufacturedLeadFieldECGVerifier::exactVm
(
    const point& p,
    const scalar t
) const
{
    switch (reference_)
    {
        case referenceType::slab:
        case referenceType::layered:
            return manufacturedLeadFieldSlabVm(p, t);
        case referenceType::sphericalShell:
            return shellTransmissionMMS::Vm(p, t, shellParams_);
    }
    return 0;
}


scalar manufacturedLeadFieldECGVerifier::exactPhi
(
    const point& p,
    const scalar t,
    const label zone
) const
{
    const leadFieldConductor& lfc = leadFieldPtr_->conductor();

    scalar phi = 0;
    scalar gauge = 0;
    switch (reference_)
    {
        case referenceType::slab:
            phi = manufacturedLeadFieldSlabPhi(p, t, slabK_);
            gauge = manufacturedLeadFieldSlabPhi(referenceCentroid_, t, slabK_);
            break;
        case referenceType::layered:
            phi = manufacturedLeadFieldLayeredPhi(p, t, slabK_, alpha_, layers_);
            gauge = manufacturedLeadFieldLayeredPhi
            (
                referenceCentroid_, t, slabK_, alpha_, layers_
            );
            break;
        case referenceType::sphericalShell:
        {
            using namespace shellTransmissionMMS;
            phi = phiE
            (
                p, t, zone < 0 ? zoneOf(p, shellParams_) : zone, shellParams_
            );
            gauge = phiE
            (
                referenceCentroid_,
                t,
                zoneOf(referenceCentroid_, shellParams_),
                shellParams_
            );
            break;
        }
    }

    // A pinned conductor fixes the constant at the pinned cell; with ground
    // patches the reference is absolute
    return
        lfc.groundPatches().empty()
      ? phi - gauge + lfc.referenceValue()
      : phi;
}


void manufacturedLeadFieldECGVerifier::solveDirect
(
    const volScalarField& Vm,
    volScalarField& phi,
    scalarField& values,
    const volScalarField* extraSource
)
{
    const leadFieldConductor& lfc = leadFieldPtr_->conductor();
    const passiveConductor& conductor = lfc.conductor();
    const dictionary& solverDict = phi.mesh().solverDict("phiE");
    const dictionary& fixedPointDict =
        phi.mesh().solutionDict().subDict("leadFieldFixedPoint");
    const scalar directTolerance = fixedPointDict.get<scalar>("tolerance");
    const label directMaxPasses = fixedPointDict.get<label>("maxPasses");
    tmp<volScalarField> tb(conductor.heartSource(Vm));
    if (extraSource)
    {
        tb.ref().primitiveFieldRef() += extraSource->primitiveField();
    }
    volScalarField& heartPhi = heartPhi_.ref();

    const scalarField& phiI = phi.primitiveField();
    scalarField phiPrev(phiI.size());
    scalar change = GREAT;
    label pass = 0;

    while (pass < 2 || change >= directTolerance)
    {
        if (pass == directMaxPasses)
        {
            FatalErrorInFunction
                << "Direct forward solve not converged after " << pass
                << " passes; last relative change " << change
                << " (tolerance " << directTolerance << ")."
                << exit(FatalError);
        }

        phiPrev = phiI;

        // The heart companion takes the base values before each assembly
        conductor.fillHeart(phi, heartPhi);

        tmp<fvScalarMatrix> tEqn(conductor.assemble(phi, heartPhi, &tb()));
        conductor.applyReference
        (
            tEqn.ref(), phi, lfc.referenceCell(), lfc.referenceValue()
        );
        solve(tEqn.ref(), solverDict);

        ++pass;
        change = gMax(mag(phiI - phiPrev))/max(gMax(mag(phiI)), VSMALL);
    }

    values.setSize(electrodeCells_.size());
    forAll(electrodeCells_, e)
    {
        values[e] = electrodeCells_[e] >= 0 ? phiI[electrodeCells_[e]] : 0;
    }
    Pstream::listCombineReduce(values, plusEqOp<scalar>());
}


void manufacturedLeadFieldECGVerifier::recordFieldError(const scalar t)
{
    const fvMesh& baseMesh = stateProvider_.baseMesh();
    const passiveConductor& conductor = leadFieldPtr_->conductor().conductor();
    const scalarField& phi = directStaticPhi_().primitiveField();
    const vectorField& C = baseMesh.C().primitiveField();
    const scalarField& V = baseMesh.V().field();

    // Heart flag with coupled neighbour values, for the interface band
    volScalarField heart
    (
        IOobject
        (
            "leadFieldVerifierHeart",
            baseMesh.time().timeName(),
            baseMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        baseMesh,
        dimensionedScalar(dimless, 0)
    );
    for (const label cellI : conductor.heartCellToBaseCell())
    {
        heart[cellI] = 1;
    }
    heart.correctBoundaryConditions();

    boolList band(baseMesh.nCells(), false);
    const labelUList& own = baseMesh.owner();
    const labelUList& nei = baseMesh.neighbour();
    forAll(nei, faceI)
    {
        if (heart[own[faceI]] != heart[nei[faceI]])
        {
            band[own[faceI]] = true;
            band[nei[faceI]] = true;
        }
    }
    forAll(heart.boundaryField(), patchI)
    {
        const fvPatchScalarField& hp = heart.boundaryField()[patchI];
        if (!hp.coupled())
        {
            continue;
        }
        const scalarField hN(hp.patchNeighbourField());
        const labelUList& cells = hp.patch().faceCells();
        forAll(cells, faceI)
        {
            if (heart[cells[faceI]] != hN[faceI])
            {
                band[cells[faceI]] = true;
            }
        }
    }

    // Cells with a face on the outer boundary of the conductor
    boolList surface(baseMesh.nCells(), false);
    forAll(baseMesh.boundary(), patchI)
    {
        const fvPatch& patch = baseMesh.boundary()[patchI];
        if (patch.coupled() || isA<emptyFvPatch>(patch))
        {
            continue;
        }
        for (const label cellI : patch.faceCells())
        {
            surface[cellI] = true;
        }
    }

    // Per region: sum |e| V, sum e^2 V, sum V, max |e|
    List<FixedList<scalar, 4>> sums(fieldRegionNames_.size(), Zero);
    scalar maxE = -1;
    point maxAt = point::zero;
    forAll(C, cellI)
    {
        const label zone =
            reference_ == referenceType::sphericalShell ? shellZones_[cellI] : -1;
        const scalar e = mag(phi[cellI] - exactPhi(C[cellI], t, zone));
        const bool inHeart = heart[cellI] > 0.5;
        const bool member[5] =
            {true, inHeart, !inHeart, band[cellI], surface[cellI]};
        for (label r = 0; r < 5; ++r)
        {
            if (member[r])
            {
                sums[r][0] += e*V[cellI];
                sums[r][1] += sqr(e)*V[cellI];
                sums[r][2] += V[cellI];
                sums[r][3] = max(sums[r][3], e);
            }
        }
        if (e > maxE)
        {
            maxE = e;
            maxAt = C[cellI];
        }
    }

    forAll(sums, r)
    {
        for (label k = 0; k < 3; ++k)
        {
            reduce(sums[r][k], sumOp<scalar>());
        }
        reduce(sums[r][3], maxOp<scalar>());
        const scalar vol = max(sums[r][2], VSMALL);
        fieldErrors_[r][0] = max(fieldErrors_[r][0], sums[r][0]/vol);
        fieldErrors_[r][1] = max(fieldErrors_[r][1], Foam::sqrt(sums[r][1]/vol));
        fieldErrors_[r][2] = max(fieldErrors_[r][2], sums[r][3]);
    }

    const scalar globalMax = returnReduce(maxE, maxOp<scalar>());
    if (maxE == globalMax)
    {
        Pout<< type() << ": largest phiStatic error " << globalMax
            << " at " << maxAt << endl;
    }
}


void manufacturedLeadFieldECGVerifier::record
(
    const List<scalar>& numericValues
)
{
    if (!leadFieldPtr_.valid())
    {
        initialise();
    }

    const volScalarField& Vm = requireVm();
    const scalar t = mesh_.time().value();
    const label nElectrodes = electrodeNames_.size();

    // Direct forward on Vm_h
    scalarField direct;
    solveDirect(Vm, directPhi_.ref(), direct);

    // Lead field and direct forward on the exact Vm at the cell centres
    volScalarField VmExact
    (
        IOobject
        (
            Vm.name(), Vm.instance(), Vm.mesh(),
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        Vm
    );
    const vectorField& centres = Vm.mesh().C().primitiveField();
    forAll(centres, cellI)
    {
        VmExact.primitiveFieldRef()[cellI] = exactVm(centres[cellI], t);
    }
    VmExact.correctBoundaryConditions();

    scalarField leadFieldStatic;
    leadFieldPtr_->evaluate(VmExact, leadFieldStatic);
    scalarField directStatic;
    if (reference_ == referenceType::sphericalShell)
    {
        const tmp<volScalarField> tS
        (
            shellSource
            (
                t,
                leadFieldPtr_->conductor().conductor()
                    .heartSource(VmExact)().dimensions()
            )
        );
        const scalarField& S = tS().primitiveField();
        const scalarField& V = stateProvider_.baseMesh().V().field();
        forAll(leadFieldStatic, e)
        {
            leadFieldStatic[e] += gSum(leadFields_[e].primitiveField()*V*S);
        }
        solveDirect(VmExact, directStaticPhi_.ref(), directStatic, &tS());
    }
    else
    {
        solveDirect(VmExact, directStaticPhi_.ref(), directStatic);
    }

    scalarField exactPoint(nElectrodes);
    scalarField exactCentroid(nElectrodes);
    forAll(exactPoint, e)
    {
        exactPoint[e] = exactPhi(electrodePositions_[e], t);
        exactCentroid[e] = exactPhi(electrodeCentroids_[e], t);
    }

    // Errors per quantity, in quantityNames_ order
    List<scalarField> errors(quantityNames_.size());
    errors[0] = scalarField(numericValues) - direct;
    errors[1] = leadFieldStatic - directStatic;
    errors[2] = leadFieldStatic - exactPoint;
    errors[3] = leadFieldStatic - exactCentroid;
    errors[4] = directStatic - exactPoint;
    errors[5] = directStatic - exactCentroid;
    errors[6] = exactPoint - exactCentroid;

    ++sampleCount_;
    forAll(errors, q)
    {
        forAll(errors[q], e)
        {
            electrodeNorms_[q][e].add(errors[q][e]);
        }
        maxNorms_[q].add(max(mag(errors[q])));
    }
    recordFieldError(t);

    relativeGapNorms_.add
    (
        max(mag(errors[0]))/max(max(mag(direct)), VSMALL)
    );
    exactScale_ = max(exactScale_, max(mag(exactPoint)));

    if (Pstream::master())
    {
        List<scalar> row;
        forAll(electrodeNames_, e)
        {
            row.append(numericValues[e]);
            row.append(direct[e]);
            row.append(exactPoint[e]);
            row.append(exactCentroid[e]);
            row.append(leadFieldStatic[e]);
            row.append(directStatic[e]);
        }
        ecgModelIO::writeRow(outputPtr_.ref(), t, row);
    }

    const Time& time = mesh_.time();
    if (time.value() + 0.5*time.deltaTValue() >= time.endTime().value())
    {
        writeSummary();
    }
}


void manufacturedLeadFieldECGVerifier::writeSummary()
{
    if (summaryWritten_ || sampleCount_ == 0)
    {
        return;
    }
    summaryWritten_ = true;

    // Reduced on every rank before the master writes
    const label nCells =
        returnReduce(stateProvider_.baseMesh().nCells(), sumOp<label>());
    const label nHeartCells = returnReduce(mesh_.nCells(), sumOp<label>());

    if (!Pstream::master())
    {
        return;
    }

    const fileName outputFile
    (
        mesh_.time().globalPath()/"postProcessing"
      / ("manufacturedLeadFieldECGSummary_" + domainName_ + ".dat")
    );
    OFstream os(outputFile);
    os.precision(std::numeric_limits<scalar>::max_digits10);

    const scalar n = sampleCount_;
    const auto writeRow = [&os, n](const word& name, const errorNorms& norms)
    {
        os  << name << " " << norms.l1/n << " " << Foam::sqrt(norms.l2/n)
            << " " << norms.linf << "\n";
    };

    os  << "# Manufactured lead-field ECG summary\n"
        << "# reference " << referenceTypeNames[reference_] << "\n";
    if (reference_ != referenceType::sphericalShell)
    {
        os  << "# k " << slabK_ << "\n";
    }
    if (reference_ == referenceType::layered)
    {
        os  << "# alpha " << alpha_ << "\n"
            << "# xLeft " << layers_.xLeft << "\n"
            << "# xRight " << layers_.xRight << "\n"
            << "# sigmaOuterLeft " << layers_.sigmaOuterLeft << "\n"
            << "# sigmaInnerLeft " << layers_.sigmaInnerLeft << "\n"
            << "# sigmaHeart " << layers_.sigmaHeart << "\n"
            << "# sigmaInnerRight " << layers_.sigmaInnerRight << "\n"
            << "# sigmaOuterRight " << layers_.sigmaOuterRight << "\n";
    }
    os  << "# cells " << nCells << "\n"
        << "# heartCells " << nHeartCells << "\n"
        << "# samples " << sampleCount_ << "\n"
        << "# time " << mesh_.time().value() << "\n"
        << "# exactScale " << exactScale_ << "\n"
        << "# Errors over the samples: <quantity> is the maximum over the "
        << "electrodes, <quantity>_<electrode> one electrode;\n"
        << "# leadField-direct_relative is that maximum divided by the "
        << "largest |direct|\n"
        << "# row L1 L2 Linf\n";

    writeRow("leadField-direct_relative", relativeGapNorms_);
    forAll(fieldRegionNames_, r)
    {
        os  << fieldRegionNames_[r] << " " << fieldErrors_[r][0] << " "
            << fieldErrors_[r][1] << " " << fieldErrors_[r][2] << "\n";
    }
    forAll(quantityNames_, q)
    {
        writeRow(quantityNames_[q], maxNorms_[q]);
        forAll(electrodeNames_, e)
        {
            writeRow
            (
                quantityNames_[q] + "_" + electrodeNames_[e],
                electrodeNorms_[q][e]
            );
        }
    }
}


void manufacturedLeadFieldECGVerifier::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
