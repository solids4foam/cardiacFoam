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

#include "leadFieldConductor.H"
#include "ecgElectrodeCells.H"
#include "myocardiumDomain.H"
#include "conductivityFieldIO.H"
#include "sealedWallPolicy.H"
#include "dimVoltage.H"
#include "HashSet.H"
#include "PstreamReduceOps.H"
#include "labelIOList.H"
#include "fvMatrices.H"
#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include <limits>

namespace Foam
{

namespace
{

word solverTypeOf(const dictionary& electroCoeffs)
{
    const word& coeffDictName = electroCoeffs.dictName();

    return coeffDictName.endsWith("Coeffs")
        ? word(coeffDictName.substr(0, coeffDictName.size() - 6))
        : coeffDictName;
}


const conductivityFieldSpec intracellularSpec
{
    "ConductivityIntracellular",
    "conductivityIntracellular",
    "conductivityIntracellular"
};

const conductivityFieldSpec extracellularSpec
{
    "ConductivityExtracellular",
    "conductivityExtracellular",
    "conductivityExtracellular"
};

const conductivityFieldSpec monodomainSpec
{
    "Conductivity",
    "conductivity",
    "conductivity"
};

// Relative tolerance of the monodomain check G_m = Gi (Gi+Ge)^-1 Ge.
const scalar monodomainConductivityTolerance = 1e-6;

// Relative tolerance on the numbers of a signature entry.
const scalar signatureTolerance = 1e-10;


// Words and labels equal; other numbers within the tolerance relative to
// the largest number of the two entries.
bool sameEntry(const tokenList& a, const tokenList& b)
{
    if (a.size() != b.size()) return false;

    scalar scale = 0.0;
    forAll(a, i)
    {
        if (a[i].isNumber()) scale = max(scale, mag(a[i].number()));
        if (b[i].isNumber()) scale = max(scale, mag(b[i].number()));
    }

    forAll(a, i)
    {
        if (a[i].isLabel() && b[i].isLabel())
        {
            if (a[i].labelToken() != b[i].labelToken()) return false;
        }
        else if (a[i].isNumber() && b[i].isNumber())
        {
            const scalar d = mag(a[i].number() - b[i].number());
            if (!(d <= signatureTolerance*scale)) return false;
        }
        else if (a[i] != b[i]) return false;
    }
    return true;
}


// Write the entries of built and current that differ to os.
bool compare(const dictionary& built, const dictionary& current, Ostream& os)
{
    wordHashSet keys(built.toc());
    keys.insert(current.toc());
    bool same = true;

    for (const word& key : keys.sortedToc())
    {
        const entry* a = built.findEntry(key, keyType::LITERAL);
        const entry* b = current.findEntry(key, keyType::LITERAL);

        if (a && b && a->isDict() && b->isDict())
        {
            same = compare(a->dict(), b->dict(), os) && same;
        }
        else if
        (
            !a || !b || a->isDict() || b->isDict()
         || !sameEntry(a->stream(), b->stream())
        )
        {
            os  << "    built: ";
            if (a) { os << *a; } else { os << key << " (none)" << nl; }
            os  << "    case:  ";
            if (b) { os << *b; } else { os << key << " (none)" << nl; }
            same = false;
        }
    }
    return same;
}


// Wall-gradient states of a potential, unless the forward's potential of
// the same name on this mesh has published them.
void publishIfAbsent(const volScalarField& psi, const bool seal)
{
    const word name(sealedWallPolicy::statesName(psi.name()));
    if (!psi.mesh().foundObject<labelIOList>(name))
    {
        sealedWallPolicy::publish(psi, seal);
    }
}


} // End anonymous namespace


leadFieldConductor::leadFieldConductor
(
    const fvMesh& baseMesh,
    const dictionary& electroCoeffs,
    const word& ecgDomainName
)
:
    baseMesh_(baseMesh),
    solverType_(solverTypeOf(electroCoeffs)),
    ownSubsetPtr_(),
    subsetPtr_(nullptr),
    GiPtr_(),
    GePtr_(),
    torsoConductivityPtr_(),
    torsoCellZones_(),
    interfaceConductivityInterpolation_(),
    sealedHeartBoundary_(electroCoeffs.get<Switch>("sealedHeartBoundary")),
    ownSealingFieldPtr_(),
    sealingFieldPtr_(nullptr),
    electrodeNames_(),
    electrodePositions_(),
    electrodeCells_(),
    groundPatches_(),
    groundValues_(),
    surfaceCurrentPatches_(),
    surfaceCurrentValues_(),
    hasReferencePoint_(false),
    referencePoint_(point::zero),
    referenceValue_(0.0),
    referenceCell_(-1),
    conductorPtr_()
{
    checkMyocardium(electroCoeffs, ecgDomainName);

    ownSubsetPtr_ = createMyocardiumMeshSubset(baseMesh_, electroCoeffs);

    if (!ownSubsetPtr_.valid() || !ownSubsetPtr_->hasSubMesh())
    {
        FatalErrorInFunction
            << "The myocardium cellZone gave no sub-mesh of mesh '"
            << baseMesh_.name() << "'."
            << exit(FatalError);
    }
    subsetPtr_ = ownSubsetPtr_.get();

    const fvMesh& heart = heartMesh();

    ownSealingFieldPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "Vm",
                heart.time().timeName(),
                heart,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            heart,
            dimensionedScalar("Vm", dimVoltage, 0.0),
            "zeroGradient"
        )
    );
    sealedWallPolicy::publish(ownSealingFieldPtr_(), sealedHeartBoundary_);
    sealingFieldPtr_ = ownSealingFieldPtr_.get();

    build(electroCoeffs, ecgDomainName);
}


leadFieldConductor::leadFieldConductor
(
    const fvMeshSubset& heartSubset,
    const volScalarField& Vm,
    const dictionary& electroCoeffs,
    const word& ecgDomainName
)
:
    baseMesh_(heartSubset.baseMesh()),
    solverType_(solverTypeOf(electroCoeffs)),
    ownSubsetPtr_(),
    subsetPtr_(&heartSubset),
    GiPtr_(),
    GePtr_(),
    torsoConductivityPtr_(),
    torsoCellZones_(),
    interfaceConductivityInterpolation_(),
    sealedHeartBoundary_(electroCoeffs.get<Switch>("sealedHeartBoundary")),
    ownSealingFieldPtr_(),
    sealingFieldPtr_(&Vm),
    electrodeNames_(),
    electrodePositions_(),
    electrodeCells_(),
    groundPatches_(),
    groundValues_(),
    surfaceCurrentPatches_(),
    surfaceCurrentValues_(),
    hasReferencePoint_(false),
    referencePoint_(point::zero),
    referenceValue_(0.0),
    referenceCell_(-1),
    conductorPtr_()
{
    checkMyocardium(electroCoeffs, ecgDomainName);

    if (&Vm.mesh() != &heartSubset.subMesh())
    {
        FatalErrorInFunction
            << "Vm is on mesh '" << Vm.mesh().name() << "', not on the "
            << "myocardium sub-mesh of mesh '" << baseMesh_.name() << "'."
            << exit(FatalError);
    }

    build(electroCoeffs, ecgDomainName);
}


void leadFieldConductor::checkMyocardium
(
    const dictionary& electroCoeffs,
    const word& ecgDomainName
) const
{
    if
    (
        solverType_ != "bidomainSolver"
     && solverType_ != "monodomainSolver"
     && solverType_ != "eikonalSolver"
    )
    {
        FatalIOErrorInFunction(electroCoeffs)
            << "ECG domain '" << ecgDomainName << "' needs a lead-field "
            << "conductor, which is defined for bidomainSolver, "
            << "monodomainSolver and eikonalSolver; myocardiumSolver is '"
            << solverType_ << "'."
            << exit(FatalIOError);
    }

    if (!electroCoeffs.found("cellZone"))
    {
        FatalIOErrorInFunction(electroCoeffs)
            << "ECG domain '" << ecgDomainName << "' needs the myocardium on "
            << "a cellZone of the heart+torso mesh; set cellZone in "
            << electroCoeffs.name() << "."
            << exit(FatalIOError);
    }
}


void leadFieldConductor::build
(
    const dictionary& electroCoeffs,
    const word& ecgDomainName
)
{
    const dictionary& ecgDict =
        electroCoeffs.subDict("ecgDomains").subDict(ecgDomainName);

    myocardiumZone_ = word(electroCoeffs.lookup("cellZone"));

    readHeartConductivities(electroCoeffs, ecgDict);

    const dictionary* bathDictPtr =
        electroCoeffs.findDict("bathPotentialDomain");
    const bool hasTorso = ecgDict.found("torso");

    if (bathDictPtr && hasTorso)
    {
        FatalIOErrorInFunction(ecgDict)
            << "ECG domain '" << ecgDomainName << "' has a torso "
            << "sub-dictionary, but " << electroCoeffs.name() << " has "
            << "bathPotentialDomain, which defines the torso conductor. "
            << "Remove one of them."
            << exit(FatalIOError);
    }

    if (!bathDictPtr && !hasTorso)
    {
        FatalIOErrorInFunction(ecgDict)
            << "ECG domain '" << ecgDomainName << "' needs a torso "
            << "sub-dictionary (bathCellZones, bathConductivityField, "
            << "interfaceConductivityInterpolation, and groundPatches or "
            << "phiERefPoint with phiEReferenceValue), because "
            << electroCoeffs.name() << " has no bathPotentialDomain."
            << exit(FatalIOError);
    }

    if (bathDictPtr)
    {
        readTorso(*bathDictPtr, true);
    }
    else
    {
        readTorso(ecgDict.subDict("torso"), false);
    }

    readElectrodes(electroCoeffs, ecgDict);
    locateElectrodes(ecgDict);

    conductorPtr_.reset
    (
        new passiveConductor
        (
            baseMesh_,
            heartMesh(),
            subsetPtr_->cellMap(),
            subsetPtr_->faceMap(),
            GiPtr_(),
            GePtr_(),
            torsoConductivityPtr_(),
            torsoCellZones_,
            interfaceConductivityInterpolation_,
            sealedHeartBoundary_,
            *sealingFieldPtr_,
            solverType_ != "bidomainSolver",
            false
        )
    );

    if (hasReferencePoint_)
    {
        referenceCell_ = conductorPtr_->referenceCell(referencePoint_);
    }
}


void leadFieldConductor::readHeartConductivities
(
    const dictionary& electroCoeffs,
    const dictionary& ecgDict
)
{
    const fvMesh& heart = heartMesh();
    const fvMeshSubset* subset = subsetPtr_;

    const wordList conductivityKeys
    ({
        "conductivitySource",
        intracellularSpec.dictionaryEntry,
        extracellularSpec.dictionaryEntry
    });

    if (solverType_ == "bidomainSolver")
    {
        forAll(conductivityKeys, keyI)
        {
            if (ecgDict.found(conductivityKeys[keyI]))
            {
                FatalIOErrorInFunction(ecgDict)
                    << "ECG domain '" << ecgDict.dictName() << "' sets "
                    << conductivityKeys[keyI] << ", but bidomainSolver "
                    << "reads Gi and Ge from " << electroCoeffs.name()
                    << ". Remove it from the ECG domain entry."
                    << exit(FatalIOError);
            }
        }

        GiPtr_.reset
        (
            readConductivityField
            (
                heart, baseMesh_, subset, electroCoeffs, intracellularSpec
            ).ptr()
        );
        GePtr_.reset
        (
            readConductivityField
            (
                heart, baseMesh_, subset, electroCoeffs, extracellularSpec
            ).ptr()
        );

        // The heart mesh may hold the myocardium's fields of these names
        GiPtr_->checkOut();
        GePtr_->checkOut();

        return;
    }

    if (!ecgDict.found("conductivitySource"))
    {
        FatalIOErrorInFunction(ecgDict)
            << "ECG domain '" << ecgDict.dictName() << "' with "
            << solverType_ << " needs Gi and Ge: conductivitySource with "
            << intracellularSpec.dictionaryEntry << " and "
            << extracellularSpec.dictionaryEntry << " (uniform), or the "
            << intracellularSpec.fieldName << " and "
            << extracellularSpec.fieldName << " fields (field)."
            << exit(FatalIOError);
    }

    GiPtr_.reset
    (
        readConductivityField
        (
            heart, baseMesh_, subset, ecgDict, intracellularSpec
        ).ptr()
    );
    GePtr_.reset
    (
        readConductivityField
        (
            heart, baseMesh_, subset, ecgDict, extracellularSpec
        ).ptr()
    );
    GiPtr_->checkOut();
    GePtr_->checkOut();

    if (solverType_ != "monodomainSolver")
    {
        return;
    }

    const tmp<volTensorField> tGm
    (
        readConductivityField
        (
            heart, baseMesh_, subset, electroCoeffs, monodomainSpec
        )
    );

    const tensorField& Gm = tGm().primitiveField();
    const tensorField& Gi = GiPtr_->primitiveField();
    const tensorField& Ge = GePtr_->primitiveField();

    forAll(Gm, cellI)
    {
        const tensor GiGe = Gi[cellI] & inv(Gi[cellI] + Ge[cellI]) & Ge[cellI];
        const scalar difference = mag(Gm[cellI] - GiGe);

        if (!(difference <= monodomainConductivityTolerance*mag(Gm[cellI])))
        {
            FatalIOErrorInFunction(ecgDict)
                << "In heart cell " << cellI << " (base cell "
                << subsetPtr_->cellMap()[cellI] << ", centre "
                << heart.C()[cellI] << ") the monodomainSolver conductivity "
                << Gm[cellI] << " differs from Gi (Gi+Ge)^-1 Ge = " << GiGe
                << " of ECG domain '" << ecgDict.dictName() << "' (Gi "
                << Gi[cellI] << ", Ge " << Ge[cellI] << "): |difference| "
                << difference << " > " << monodomainConductivityTolerance
                << " |conductivity|."
                << exit(FatalIOError);
        }
    }
}


void leadFieldConductor::readTorso
(
    const dictionary& torsoDict,
    const bool bathBlock
)
{
    if (!bathBlock)
    {
        const wordHashSet validKeys
        ({
            "bathCellZones",
            "bathConductivityField",
            "interfaceConductivityInterpolation",
            "groundPatches",
            "surfaceCurrentPatches",
            "phiERefPoint",
            "phiEReferenceValue"
        });

        const wordList keys(torsoDict.toc());

        forAll(keys, keyI)
        {
            if (!validKeys.found(keys[keyI]))
            {
                FatalIOErrorInFunction(torsoDict)
                    << "Unknown key '" << keys[keyI] << "' in "
                    << torsoDict.name() << ". Valid keys are "
                    << validKeys.sortedToc() << "."
                    << exit(FatalIOError);
            }
        }
    }

    const passiveConductor::torsoSettings torso =
        passiveConductor::torsoSettings::read(baseMesh_, torsoDict);

    torsoCellZones_ = torso.bathCellZones;
    interfaceConductivityInterpolation_ =
        torso.interfaceConductivityInterpolation;
    torsoConductivityPtr_ = passiveConductor::readTorsoConductivity
    (
        baseMesh_, torso.bathConductivityField, false
    );
    groundPatches_ = torso.groundPatches;
    groundValues_ = torso.groundValues;
    surfaceCurrentPatches_ = torso.surfaceCurrentPatches;
    surfaceCurrentValues_ = torso.surfaceCurrentValues;

    // Ground patches fix the constant; a reference point is then ignored,
    // as in the forward. Otherwise the constant is pinned.
    if (groundPatches_.empty())
    {
        torso.requirePhiERefPoint("leadFieldConductor");
        hasReferencePoint_ = true;
        referencePoint_ = torso.phiERefPoint;
        referenceValue_ = torso.phiEReferenceValue;
    }
}


void leadFieldConductor::readElectrodes
(
    const dictionary& electroCoeffs,
    const dictionary& ecgDict
)
{
    const dictionary* electrodeDictPtr = ecgDict.findDict("electrodePositions");

    if (!electrodeDictPtr)
    {
        electrodeDictPtr =
            electroCoeffs.subDict("ecgDomains").findDict("electrodePositions");
    }

    if (!electrodeDictPtr || electrodeDictPtr->empty())
    {
        FatalIOErrorInFunction(ecgDict)
            << "ECG domain '" << ecgDict.dictName() << "' has no electrodes; "
            << "set electrodePositions in the ECG domain entry or in "
            << "ecgDomains."
            << exit(FatalIOError);
    }

    electrodeNames_ = electrodeDictPtr->toc();
    electrodePositions_.setSize(electrodeNames_.size());

    forAll(electrodeNames_, electrodeI)
    {
        electrodePositions_[electrodeI] =
            electrodeDictPtr->get<vector>(electrodeNames_[electrodeI]);
    }
}


void leadFieldConductor::locateElectrodes(const dictionary& ecgDict)
{
    electrodeCells_ = ecgElectrodeCells(baseMesh_, electrodePositions_);

    boolList inHeart(baseMesh_.nCells(), false);
    for (const label cellI : subsetPtr_->cellMap())
    {
        inHeart[cellI] = true;
    }

    forAll(electrodeCells_, electrodeI)
    {
        const label cellI = electrodeCells_[electrodeI];

        if (returnReduce(cellI >= 0 && inHeart[cellI], orOp<bool>()))
        {
            FatalIOErrorInFunction(ecgDict)
                << "Electrode '" << electrodeNames_[electrodeI] << "' at "
                << electrodePositions_[electrodeI] << " of ECG domain '"
                << ecgDict.dictName() << "' lies in a cell of the myocardium "
                << "cellZone '" << myocardiumZone_ << "'. Lead-field "
                << "electrodes must lie in the torso."
                << exit(FatalIOError);
        }
    }
}


dictionary leadFieldConductor::signature() const
{
    // Entries keep every digit, so that sums compare at a relative 1e-10
    const unsigned precision = IOstream::defaultPrecision
    (
        std::numeric_limits<scalar>::max_digits10
    );

    dictionary sig;

    sig.add("nBaseCells", returnReduce(baseMesh_.nCells(), sumOp<label>()));
    sig.add("nHeartCells", returnReduce(heartMesh().nCells(), sumOp<label>()));
    sig.add("sumBaseVolume", gSum(baseMesh_.V().field()));
    sig.add
    (
        "sumBaseVolumeCentre",
        gSum
        (
            vectorField
            (
                baseMesh_.V().field()*baseMesh_.C().primitiveField()
            )
        )
    );

    // Zone ids looked up once: the myocardium zone, then the torso zones.
    // The passiveConductor constructor has already made a missing torso
    // zone fatal (the single validation point), so only the myocardium
    // zone can be absent here.
    wordList zoneNames(1, myocardiumZone_);
    zoneNames.append(torsoCellZones_);

    labelList zoneIds(zoneNames.size());
    forAll(zoneNames, zoneI)
    {
        zoneIds[zoneI] = baseMesh_.cellZones().findZoneID(zoneNames[zoneI]);
    }

    dictionary zoneCounts;
    {
        forAll(zoneNames, zoneI)
        {
            const label zoneId = zoneIds[zoneI];

            const label nZoneCells =
                zoneId < 0 ? 0 : baseMesh_.cellZones()[zoneId].size();

            zoneCounts.add
            (
                zoneNames[zoneI],
                returnReduce(nZoneCells, sumOp<label>())
            );
        }
    }
    sig.add("zoneCellCounts", zoneCounts);

    sig.add("sumGi", gSum(GiPtr_->primitiveField()));
    sig.add("sumGe", gSum(GePtr_->primitiveField()));

    boolList inTorso(baseMesh_.nCells(), false);
    forAll(torsoCellZones_, zoneI)
    {
        // Torso zone zoneI sits after the myocardium zone in zoneIds
        for (const label cellI : baseMesh_.cellZones()[zoneIds[zoneI + 1]])
        {
            inTorso[cellI] = true;
        }
    }

    const scalarField& torsoSigma = torsoConductivityPtr_->primitiveField();
    scalar sumTorso = 0.0;
    forAll(inTorso, cellI)
    {
        if (inTorso[cellI])
        {
            sumTorso += torsoSigma[cellI];
        }
    }
    reduce(sumTorso, sumOp<scalar>());
    sig.add("sumTorsoConductivity", sumTorso);

    sig.add
    (
        "interfaceConductivityInterpolation",
        interfaceConductivityInterpolation_
    );
    sig.add("sealedHeartBoundary", sealedHeartBoundary_);

    dictionary sealingTypes;
    const volScalarField::Boundary& sealingBf =
        sealingFieldPtr_->boundaryField();
    forAll(sealingBf, patchI)
    {
        if (!sealingBf[patchI].coupled())
        {
            sealingTypes.add
            (
                sealingBf[patchI].patch().name(),
                sealingBf[patchI].type()
            );
        }
    }
    sig.add("sealingPatchTypes", sealingTypes);

    dictionary electrodes;
    forAll(electrodeNames_, electrodeI)
    {
        electrodes.add
        (
            electrodeNames_[electrodeI],
            electrodePositions_[electrodeI]
        );
    }
    sig.add("electrodePositions", electrodes);

    dictionary ground;
    forAll(groundPatches_, patchI)
    {
        ground.add(groundPatches_[patchI], groundValues_[patchI]);
    }
    sig.add("groundPatches", ground);

    dictionary currents;
    forAll(surfaceCurrentPatches_, patchI)
    {
        currents.add
        (
            surfaceCurrentPatches_[patchI],
            surfaceCurrentValues_[patchI]
        );
    }
    sig.add("surfaceCurrentPatches", currents);

    if (hasReferencePoint_)
    {
        sig.add("phiERefPoint", referencePoint_);
    }
    sig.add("phiEReferenceValue", referenceValue_);

    IOstream::defaultPrecision(precision);

    return sig;
}

bool leadFieldConductor::matches(const dictionary& built, Ostream& os) const
{
    return compare(built, signature(), os);
}


tmp<volScalarField> leadFieldConductor::forwardPotential() const
{
    const fvMesh& mesh = baseMesh_;

    wordList types
    (
        mesh.boundary().size(), zeroGradientFvPatchScalarField::typeName
    );
    for (const word& patch : groundPatches_)
    {
        types[mesh.boundaryMesh().findPatchID(patch)] =
            fixedValueFvPatchScalarField::typeName;
    }
    for (const word& patch : surfaceCurrentPatches_)
    {
        types[mesh.boundaryMesh().findPatchID(patch)] =
            fixedGradientFvPatchScalarField::typeName;
    }

    tmp<volScalarField> tphiE
    (
        new volScalarField
        (
            IOobject
            (
                "phiE", mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimVoltage, Zero),
            types
        )
    );
    volScalarField& phiE = tphiE.ref();

    forAll(groundPatches_, i)
    {
        const label patchI = mesh.boundaryMesh().findPatchID(groundPatches_[i]);
        phiE.boundaryFieldRef()[patchI] == groundValues_[i];
    }

    // Surface current I as the forward imposes it: gradient I/(n.sigma.n)
    forAll(surfaceCurrentPatches_, i)
    {
        const label patchI =
            mesh.boundaryMesh().findPatchID(surfaceCurrentPatches_[i]);
        auto& pf = refCast<fixedGradientFvPatchScalarField>
        (
            phiE.boundaryFieldRef()[patchI]
        );
        const tensorField sigma
        (
            conductorPtr_->sigmaTotal().boundaryField()[patchI]
           .patchInternalField()
        );
        const vectorField n(pf.patch().nf());
        pf.gradient() = surfaceCurrentValues_[i]/(n & sigma & n);
        pf.evaluate();
    }

    publishIfAbsent(phiE, false);

    return tphiE;
}


tmp<volScalarField> leadFieldConductor::heartPotential() const
{
    const fvMesh& heart = heartMesh();

    tmp<volScalarField> theartPhiE
    (
        new volScalarField
        (
            IOobject
            (
                "phiE", heart.time().timeName(), heart,
                IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
            ),
            heart,
            dimensionedScalar(dimVoltage, Zero),
            "zeroGradient"
        )
    );
    publishIfAbsent(theartPhiE(), sealedHeartBoundary_);

    return theartPhiE;
}


tmp<scalarField> leadFieldConductor::constantSource() const
{
    const passiveConductor& conductor = conductorPtr_();
    tmp<volScalarField> tphiE(forwardPotential());
    tmp<volScalarField> theartPhiE(heartPotential());
    volScalarField& phiE = tphiE.ref();

    tmp<fvScalarMatrix> tEqn(conductor.assemble(phiE, theartPhiE.ref()));
    conductor.applyReference(tEqn.ref(), phiE, referenceCell_, referenceValue_);

    // The source the linear solver sees: the matrix source plus the
    // non-coupled boundary coefficients
    tmp<scalarField> tb(new scalarField(tEqn().source()));
    forAll(phiE.boundaryField(), patchI)
    {
        if (phiE.boundaryField()[patchI].coupled()) continue;

        const labelUList& cells = baseMesh_.boundary()[patchI].faceCells();
        const scalarField& coeffs = tEqn().boundaryCoeffs()[patchI];
        forAll(cells, faceI)
        {
            tb.ref()[cells[faceI]] += coeffs[faceI];
        }
    }

    return tb;
}

} // End namespace Foam

// ************************************************************************* //
