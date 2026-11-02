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
#include "torsoCoupling.H"
#include "dimVoltage.H"
#include "HashSet.H"
#include "PstreamReduceOps.H"
#include "labelIOList.H"
#include "fvMatrices.H"
#include "fvGeometryScheme.H"
#include "zeroGradientFvPatchFields.H"
#include <limits>

namespace Foam
{

namespace
{

// The <myocardiumSolver>Coeffs block of electroProperties
const dictionary& myocardiumCoeffs(const dictionary& electroProperties)
{
    return electroProperties.subDict
    (
        electroProperties.get<word>("myocardiumSolver") + "Coeffs"
    );
}


// The entry ecgDomainName of the ecgDomains block of electroProperties
const dictionary& ecgDomainDict
(
    const dictionary& electroProperties,
    const word& ecgDomainName
)
{
    return electroProperties.subDict("ecgDomains").subDict(ecgDomainName);
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


// The fvSchemes entries that can resolve a scheme of phiE: default,
// patterns and keys containing phiE, of the laplacian, snGrad, gradient and
// interpolation schemes. Keys are written as words.
dictionary phiESchemes(const fvMesh& mesh)
{
    dictionary schemes;
    for
    (
        const word& group
      : {
            "laplacianSchemes", "snGradSchemes", "gradSchemes",
            "interpolationSchemes"
        }
    )
    {
        dictionary entries;
        const dictionary* groupDict = mesh.schemesDict().findDict(group);
        if (groupDict)
        {
            for (const entry& e : *groupDict)
            {
                const keyType& key = e.keyword();
                if
                (
                    e.isStream()
                 && (
                        key == "default" || key.isPattern()
                     || key.find("phiE") != std::string::npos
                    )
                )
                {
                    entries.add(word::validate(key), e.stream());
                }
            }
        }
        schemes.add(group, entries);
    }
    return schemes;
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


leadFieldConductor::fixedPointControls
leadFieldConductor::fixedPointControls::read(const fvMesh& mesh)
{
    const dictionary& dict =
        mesh.solutionDict().subDict("leadFieldFixedPoint");

    fixedPointControls controls;
    controls.tolerance = dict.get<scalar>("tolerance");
    controls.maxPasses = dict.get<label>("maxPasses");
    return controls;
}


leadFieldConductor::leadFieldConductor
(
    const fvMesh& baseMesh,
    const dictionary& electroProperties,
    const word& ecgDomainName,
    const ecgElectrodes& electrodes
)
:
    baseMesh_(baseMesh),
    solverType_(electroProperties.get<word>("myocardiumSolver")),
    bidomain_(solverType_ == "bidomainSolver"),
    bathHeartPhiETrace_(bidomainSolver::heartPhiETrace::zeroGradient),
    sealedHeartBoundary_
    (
        myocardiumCoeffs(electroProperties).get<Switch>("sealedHeartBoundary")
    ),
    myocardiumZone_(),
    electrodes_(electrodes),
    ownSubsetPtr_(),
    subsetPtr_(nullptr),
    ownGiPtr_(),
    ownGePtr_(),
    GiPtr_(nullptr),
    GePtr_(nullptr),
    ownSealingFieldPtr_(),
    sealingFieldPtr_(nullptr),
    electrodeCells_(),
    conductorPtr_()
{
    checkMyocardium(electroProperties, ecgDomainName);

    const dictionary& electroCoeffs = myocardiumCoeffs(electroProperties);
    myocardiumZone_ = electroCoeffs.get<word>("cellZone");
    if (bidomain_)
    {
        bathHeartPhiETrace_ =
            bidomainSolver::readBathHeartPhiETrace(electroCoeffs);

        if (bathHeartPhiETrace_ != bidomainSolver::heartPhiETrace::zeroGradient)
        {
            FatalIOErrorInFunction(electroCoeffs)
                << "ECG domain '" << ecgDomainName << "': lead fields are "
                << "built with the zeroGradient heart phiE companion; "
                << "bathHeartPhiETrace "
                << bidomainSolver::heartPhiETraceNames[bathHeartPhiETrace_]
                << " is not supported."
                << exit(FatalIOError);
        }
    }
}


leadFieldConductor::leadFieldConductor
(
    const fvMesh& baseMesh,
    const dictionary& electroProperties,
    const word& ecgDomainName
)
:
    leadFieldConductor
    (
        baseMesh,
        electroProperties,
        ecgDomainName,
        readEcgElectrodes
        (
            ecgDomainDict(electroProperties, ecgDomainName),
            electroProperties.subDict("ecgDomains")
           .findDict("electrodePositions")
        )
    )
{
    ownSubsetPtr_ =
        createMyocardiumMeshSubset
        (
            baseMesh_, myocardiumCoeffs(electroProperties)
        );

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

    build(electroProperties, ecgDomainName, nullptr);
}


leadFieldConductor::leadFieldConductor
(
    const electroStateProvider& heart,
    const dictionary& electroProperties,
    const word& ecgDomainName,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    leadFieldConductor
    (
        heart.baseMesh(),
        electroProperties,
        ecgDomainName,
        ecgElectrodes{electrodeNames, electrodePositions}
    )
{
    subsetPtr_ = heart.meshSubsetPtr();
    sealingFieldPtr_ = heart.VmPtr();

    if (!subsetPtr_ || !subsetPtr_->hasSubMesh() || !sealingFieldPtr_)
    {
        FatalErrorInFunction
            << "ECG domain '" << ecgDomainName << "' needs the myocardium "
            << "sub-mesh of its cellZone '" << myocardiumZone_ << "' and its "
            << "Vm."
            << exit(FatalError);
    }

    if (&sealingFieldPtr_->mesh() != &subsetPtr_->subMesh())
    {
        FatalErrorInFunction
            << "Vm is on mesh '" << sealingFieldPtr_->mesh().name()
            << "', not on the myocardium sub-mesh of mesh '"
            << baseMesh_.name() << "'."
            << exit(FatalError);
    }

    build(electroProperties, ecgDomainName, &heart);
}


void leadFieldConductor::checkMyocardium
(
    const dictionary& electroProperties,
    const word& ecgDomainName
) const
{
    const torsoCoupling::type coupling = torsoCoupling::read(electroProperties);

    if (coupling != torsoCoupling::type::leadField)
    {
        FatalIOErrorInFunction(electroProperties)
            << "ECG domain '" << ecgDomainName << "' needs lead fields, "
            << "which are built for torsoCoupling leadField; torsoCoupling is "
            << torsoCoupling::name(coupling) << "."
            << (
                   coupling == torsoCoupling::type::bath
                 ? " With torsoCoupling bath, torsoECG samples the phiE that "
                   "the bath solves and needs no lead fields."
                 : ""
               )
            << exit(FatalIOError);
    }

    const dictionary& electroCoeffs = myocardiumCoeffs(electroProperties);

    if (!bidomain_ && solverType_ != "monodomainSolver")
    {
        FatalIOErrorInFunction(electroCoeffs)
            << "ECG domain '" << ecgDomainName << "' needs a lead-field "
            << "conductor, which is defined for bidomainSolver and "
            << "monodomainSolver; myocardiumSolver is '" << solverType_
            << "'. eikonalSolver is not supported."
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
    const dictionary& electroProperties,
    const word& ecgDomainName,
    const electroStateProvider* heart
)
{
    const dictionary& ecgDict =
        ecgDomainDict(electroProperties, ecgDomainName);

    setHeartConductivities
    (
        myocardiumCoeffs(electroProperties), ecgDict, heart
    );
    locateElectrodes(ecgDict);

    // A monodomain source is sealed at the heart wall, so the heart side of
    // heart-torso faces carries Gi + Ge; a bidomain heart carries Ge
    const bool totalWallConductivity = !bidomain_;

    // The conductor's base-mesh fields stay out of the base-mesh registry
    const bool registerFields = false;

    conductorPtr_.reset
    (
        new passiveConductor
        (
            baseMesh_,
            heartMesh(),
            subsetPtr_->cellMap(),
            subsetPtr_->faceMap(),
            *GiPtr_,
            *GePtr_,
            passiveConductor::torsoSettings::read
            (
                baseMesh_, electroProperties.subDict("torso")
            ),
            sealedHeartBoundary_,
            *sealingFieldPtr_,
            totalWallConductivity,
            registerFields
        )
    );
}


void leadFieldConductor::setHeartConductivities
(
    const dictionary& electroCoeffs,
    const dictionary& ecgDict,
    const electroStateProvider* heart
)
{
    const fvMesh& heartM = heartMesh();
    const fvMeshSubset* subset = subsetPtr_;

    if (bidomain_)
    {
        for
        (
            const word& key
          : {
                word("conductivitySource"),
                intracellularSpec.dictionaryEntry,
                extracellularSpec.dictionaryEntry
            }
        )
        {
            if (ecgDict.found(key))
            {
                FatalIOErrorInFunction(ecgDict)
                    << "ECG domain '" << ecgDict.dictName() << "' sets "
                    << key << ", but bidomainSolver reads Gi and Ge from "
                    << electroCoeffs.name() << ". Remove it from the ECG "
                    << "domain entry."
                    << exit(FatalIOError);
            }
        }

        if (heart)
        {
            // The conductivities of the running bidomain myocardium
            GiPtr_ = heart->intracellularConductivityPtr();
            GePtr_ = heart->extracellularConductivityPtr();

            if
            (
                !GiPtr_ || !GePtr_
             || &GiPtr_->mesh() != &heartM || &GePtr_->mesh() != &heartM
            )
            {
                FatalErrorInFunction
                    << "The bidomain myocardium exposes no intracellular and "
                    << "extracellular conductivities on its sub-mesh."
                    << exit(FatalError);
            }
            return;
        }

        ownGiPtr_.reset
        (
            readConductivityField
            (
                heartM, baseMesh_, subset, electroCoeffs, intracellularSpec
            ).ptr()
        );
        ownGePtr_.reset
        (
            readConductivityField
            (
                heartM, baseMesh_, subset, electroCoeffs, extracellularSpec
            ).ptr()
        );
    }
    else
    {
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

        ownGiPtr_.reset
        (
            readConductivityField
            (
                heartM, baseMesh_, subset, ecgDict, intracellularSpec
            ).ptr()
        );
        ownGePtr_.reset
        (
            readConductivityField
            (
                heartM, baseMesh_, subset, ecgDict, extracellularSpec
            ).ptr()
        );
    }

    // The heart mesh may hold the myocardium's fields of these names
    ownGiPtr_->checkOut();
    ownGePtr_->checkOut();
    GiPtr_ = ownGiPtr_.get();
    GePtr_ = ownGePtr_.get();

    if (bidomain_)
    {
        return;
    }

    const tmp<volTensorField> tGm
    (
        readConductivityField
        (
            heartM, baseMesh_, subset, electroCoeffs, monodomainSpec
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
                << heartM.C()[cellI] << ") the monodomainSolver conductivity "
                << Gm[cellI] << " differs from Gi (Gi+Ge)^-1 Ge = " << GiGe
                << " of ECG domain '" << ecgDict.dictName() << "' (Gi "
                << Gi[cellI] << ", Ge " << Ge[cellI] << "): |difference| "
                << difference << " > " << monodomainConductivityTolerance
                << " |conductivity|."
                << exit(FatalIOError);
        }
    }
}


void leadFieldConductor::locateElectrodes(const dictionary& ecgDict)
{
    electrodeCells_ =
        ecgElectrodeCells(baseMesh_, electrodes_.names, electrodes_.positions);

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
                << "Electrode '" << electrodes_.names[electrodeI] << "' at "
                << electrodes_.positions[electrodeI] << " of ECG domain '"
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

    sig.add("operator", word("adjoint"));
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

    const passiveConductor& conductor = conductorPtr_();
    const passiveConductor::torsoSettings& torso = conductor.torso();

    // The myocardium zone, then the torso zones; the torso zones exist
    // (checked by passiveConductor), the myocardium zone may not
    wordList zoneNames(1, myocardiumZone_);
    zoneNames.append(torso.cellZones);

    dictionary zoneCounts;
    for (const word& zoneName : zoneNames)
    {
        const label zoneId = baseMesh_.cellZones().findZoneID(zoneName);
        const label nZoneCells =
            zoneId < 0 ? 0 : baseMesh_.cellZones()[zoneId].size();

        zoneCounts.add(zoneName, returnReduce(nZoneCells, sumOp<label>()));
    }
    sig.add("zoneCellCounts", zoneCounts);

    sig.add("sumGi", gSum(GiPtr_->primitiveField()));
    sig.add("sumGe", gSum(GePtr_->primitiveField()));

    // The torso cells are the cells outside the heart (every cell is in the
    // myocardium or a torso zone, checked by passiveConductor), where
    // sigmaTotal is the torso conductivity times the identity
    boolList inHeart(baseMesh_.nCells(), false);
    for (const label cellI : conductor.heartCellToBaseCell())
    {
        inHeart[cellI] = true;
    }

    const tensorField& sigmaTotal = conductor.sigmaTotal().primitiveField();
    scalar sumTorso = 0.0;
    forAll(inHeart, cellI)
    {
        if (!inHeart[cellI])
        {
            sumTorso += sigmaTotal[cellI].xx();
        }
    }
    reduce(sumTorso, sumOp<scalar>());
    sig.add("sumTorsoConductivity", sumTorso);

    sig.add("myocardiumSolver", solverType_);
    if (bidomain_)
    {
        sig.add
        (
            "bathHeartPhiETrace",
            bidomainSolver::heartPhiETraceNames[bathHeartPhiETrace_]
        );
    }
    sig.add("geometryScheme", baseMesh_.geometry().type());
    sig.add("baseSchemes", phiESchemes(baseMesh_));
    sig.add("heartSchemes", phiESchemes(heartMesh()));
    sig.add
    (
        "interfaceConductivityInterpolation",
        passiveConductor::interfaceInterpolationNames[torso.interface]
    );
    if
    (
        torso.interface
     == passiveConductor::interfaceInterpolation::conormalHarmonic
    )
    {
        sig.add
        (
            "nInterfaceGradientCorrectors",
            conductor.nInterfaceGradientCorrectors()
        );
    }
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
    forAll(electrodes_.names, electrodeI)
    {
        electrodes.add
        (
            electrodes_.names[electrodeI],
            electrodes_.positions[electrodeI]
        );
    }
    sig.add("electrodePositions", electrodes);

    dictionary ground;
    forAll(torso.groundPatches, patchI)
    {
        ground.add(torso.groundPatches[patchI], torso.groundValues[patchI]);
    }
    sig.add("groundPatches", ground);

    dictionary currents;
    forAll(torso.surfaceCurrentPatches, patchI)
    {
        currents.add
        (
            torso.surfaceCurrentPatches[patchI],
            torso.surfaceCurrentValues[patchI]
        );
    }
    sig.add("surfaceCurrentPatches", currents);

    // The reference of a pinned conductor; ground patches fix the constant
    if (torso.pinsReference())
    {
        sig.add("phiERefPoint", torso.phiERefPoint);
    }
    sig.add
    (
        "phiEReferenceValue",
        torso.pinsReference() ? torso.phiEReferenceValue : 0.0
    );

    IOstream::defaultPrecision(precision);

    return sig;
}

bool leadFieldConductor::matches(const dictionary& built, Ostream& os) const
{
    return compare(built, signature(), os);
}


tmp<volScalarField> leadFieldConductor::newPotential
(
    const fvMesh& mesh,
    const dimensionSet& dims,
    const wordList& patchTypes,
    const bool seal
) const
{
    tmp<volScalarField> tpsi
    (
        new volScalarField
        (
            IOobject
            (
                passiveConductor::potentialName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dims, Zero),
            patchTypes
        )
    );
    publishIfAbsent(tpsi(), seal);

    return tpsi;
}


tmp<volScalarField> leadFieldConductor::forwardPotential() const
{
    const passiveConductor& conductor = conductorPtr_();

    tmp<volScalarField> tphiE
    (
        newPotential
        (
            baseMesh_,
            dimVoltage,
            conductor.torso().patchTypes(baseMesh_, true),
            false
        )
    );
    conductor.torso().setGroundValues(tphiE.ref());
    conductor.setSurfaceCurrents(tphiE.ref());

    return tphiE;
}


tmp<volScalarField> leadFieldConductor::leadFieldPotential() const
{
    const passiveConductor& conductor = conductorPtr_();

    // z is in 1/(conductivity length): A z = e in units of the source
    return newPotential
    (
        baseMesh_,
        dimless/(conductor.sigmaTotal().dimensions()*dimLength),
        conductor.torso().patchTypes(baseMesh_, false),
        false
    );
}


tmp<volScalarField> leadFieldConductor::heartPotential() const
{
    const fvMesh& heart = heartMesh();

    return newPotential
    (
        heart,
        dimVoltage,
        wordList
        (
            heart.boundary().size(), zeroGradientFvPatchScalarField::typeName
        ),
        sealedHeartBoundary_
    );
}


tmp<scalarField> leadFieldConductor::constantSource() const
{
    const passiveConductor& conductor = conductorPtr_();
    tmp<volScalarField> tphiE(forwardPotential());
    tmp<volScalarField> theartPhiE(heartPotential());
    volScalarField& phiE = tphiE.ref();

    tmp<fvScalarMatrix> tEqn(conductor.assemble(phiE, theartPhiE.ref()));
    conductor.applyReference
    (
        tEqn.ref(), conductor.torso().phiEReferenceValue
    );

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
