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

#include "leadFieldData.H"
#include "conductivityFieldIO.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "cellZoneMesh.H"
#include "sealedWallPolicy.H"
#include "labelIOList.H"
#include "gradScheme.H"
#include "PstreamReduceOps.H"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace Foam
{
namespace
{
const dimensionSet sigmaDims(-1, -3, 3, 0, 0, 2, 0);
const dimensionSet sensitivityDims(0, -2, 0, 0, 0, 0, 0);

void validatePhysical(const volTensorField& f)
{
    if (f.dimensions() != sigmaDims)
    {
        FatalErrorInFunction << f.name() << " must have conductivity dimensions S/m."
            << exit(FatalError);
    }
    bool valid = true;
    for (const tensor& g : f.primitiveField())
    {
        for (direction d = 0; d < tensor::nComponents; ++d)
            valid = valid && std::isfinite(g[d]);
        const scalar scale = max(mag(g), VSMALL);
        const scalar tol = 1e-12*scale;
        valid = valid && mag(g-g.T()) <= tol
            && g.xx() >= -tol && g.yy() >= -tol && g.zz() >= -tol
            && g.xx()*g.yy()-sqr(g.xy()) >= -tol*scale
            && g.xx()*g.zz()-sqr(g.xz()) >= -tol*scale
            && g.yy()*g.zz()-sqr(g.yz()) >= -tol*scale
            && det(g) >= -tol*sqr(scale);
    }
    reduce(valid, andOp<bool>());
    if (!valid)
    {
        FatalErrorInFunction << f.name()
            << " must contain finite symmetric positive-semidefinite conductivity tensors."
            << exit(FatalError);
    }
}

template<class FieldType>
tmp<FieldType> readStatic
(
    const fvMesh& mesh, const word& name, const dictionary& dict,
    const dimensionSet& dims
)
{
    tmp<FieldType> t
    (
        new FieldType
        (
            IOobject(name, dict.lookupOrDefault<word>("fieldInstance", "0"),
                mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
        )
    );
    if (t().dimensions() != dims)
    {
        FatalErrorInFunction << "Lead-field input '" << name
            << "' has wrong dimensions; expected " << dims << exit(FatalError);
    }
    return t;
}

label sourceCell(const electroStateProvider& state, label cell)
{
    const labelUList* map = state.subsetCellMapPtr();
    if (&state.mesh() == &state.baseMesh()) return cell;
    if (!map || map->size() != state.mesh().nCells()
     || (*map)[cell] < 0 || (*map)[cell] >= state.baseMesh().nCells())
    {
        FatalErrorInFunction << "Lead-field import requires a valid myocardial "
            << "subset-to-base cell map." << exit(FatalError);
    }
    return (*map)[cell];
}

bool sameDictionary(const dictionary& a, const dictionary& b)
{
    if (a.size()!=b.size()) return false;
    forAllConstIter(dictionary,a,iter)
    {
        const entry& x=iter();
        if (!b.found(x.keyword())) return false;
        const entry& y=b.lookupEntry(x.keyword(),keyType::LITERAL);
        if (x.isDict()!=y.isDict()) return false;
        if (x.isDict())
        {
            if (!sameDictionary(x.dict(),y.dict())) return false;
        }
        else
        {
            const ITstream& u=x.stream();
            const ITstream& v=y.stream();
            if (u.size()!=v.size()) return false;
            forAll(u,i)
            {
                if (u[i].isNumber() && v[i].isNumber())
                {
                    if (u[i].number()!=v[i].number()) return false;
                }
                else if (u[i]!=v[i]) return false;
            }
        }
    }
    return true;
}

void hashValue(std::uint64_t& h, std::uint64_t value)
{
    for (unsigned i=0; i<8; ++i)
    {
        h ^= (value >> (8*i)) & 255;
        h *= UINT64_C(1099511628211);
    }
}
}

tmp<volTensorField> leadFieldData::readPhysical
(
    const fvMesh& mesh, const dictionary& dict, bool intracellular
)
{
    // A legacy provider pointer is insufficient: require an explicit source.
    (void)dict.get<word>("conductivitySource");
    const conductivityFieldSpec spec = intracellular
      ? conductivityFieldSpec{"ConductivityIntracellular", "conductivityIntracellular", "leadFieldGiInput", "conductivityIntracellular"}
      : conductivityFieldSpec{"ConductivityExtracellular", "conductivityExtracellular", "leadFieldGeInput", "conductivityExtracellular"};
    tmp<volTensorField> t = readConductivityField(mesh, mesh, nullptr, dict, spec);
    validatePhysical(t());
    return t;
}

vectorField leadFieldData::meshSignatures(const fvMesh& mesh)
{
    vectorField result(mesh.nCells());
    forAll(result, cellI)
    {
        std::vector<std::uint64_t> faceHashes;
        for (label faceI : mesh.cells()[cellI])
        {
            std::vector<point> pts;
            for (label pointI : mesh.faces()[faceI]) pts.push_back(mesh.points()[pointI]);
            std::sort(pts.begin(), pts.end(), [](const point& a, const point& b)
            {
                if (a.x()!=b.x()) return a.x()<b.x();
                if (a.y()!=b.y()) return a.y()<b.y();
                return a.z()<b.z();
            });
            std::uint64_t h = UINT64_C(14695981039346656037);
            hashValue(h, pts.size());
            for (const point& p : pts) for (direction d=0; d<3; ++d)
            {
                const double value = p[d] == 0 ? 0 : p[d];
                std::uint64_t bits;
                std::memcpy(&bits, &value, sizeof(bits));
                hashValue(h, bits);
            }
            faceHashes.push_back(h);
        }
        std::sort(faceHashes.begin(), faceHashes.end());
        std::uint64_t h = UINT64_C(14695981039346656037);
        for (auto f : faceHashes) hashValue(h, f);
        result[cellI] = vector(scalar(h & UINT64_C(0xffff)), scalar((h >> 16) & UINT64_C(0xffff)), scalar((h >> 32) & UINT64_C(0xffff)));
    }
    return result;
}

void leadFieldData::writeProvenance
(
    const fvMesh& mesh, const volTensorField& Gi, const volTensorField& Ge,
    const dictionary& props, const word& gradientConstruction
)
{
    const word instance(props.lookupOrDefault<word>("fieldInstance", "0"));
    volTensorField gi(IOobject("leadFieldGi", instance, mesh, IOobject::NO_READ,
        IOobject::NO_WRITE), Gi);
    volTensorField ge(IOobject("leadFieldGe", instance, mesh, IOobject::NO_READ,
        IOobject::NO_WRITE), Ge);
    volVectorField signature(IOobject("leadFieldMeshSignature", instance, mesh,
        IOobject::NO_READ, IOobject::NO_WRITE), mesh,
        dimensionedVector("zero", dimless, Zero));
    signature.primitiveFieldRef() = meshSignatures(mesh);
    gi.write(); ge.write(); signature.write();
    IOdictionary metadata(IOobject("leadFieldMetadata", mesh.time().constant(),
        mesh, IOobject::NO_READ, IOobject::NO_WRITE));
    metadata.add("formatVersion", 1);
    metadata.add("representation", word("conductivityWeightedVector"));
    metadata.add("gradientConstruction",gradientConstruction);
    metadata.add("leadFields", props.subDict("leadFields"));
    metadata.add("inputProperties", props);
    metadata.add("supportCells", returnReduce(mesh.nCells(),sumOp<label>()));
    metadata.add("geometricDimensions", mesh.nGeometricD());
    label heartCells=mesh.nCells();
    if (props.found("cellZone"))
    {
        const label zone=mesh.cellZones().findZoneID(props.get<word>("cellZone"));
        if (zone<0) FatalErrorInFunction << "Unknown lead-field heart cellZone." << exit(FatalError);
        volScalarField region(IOobject("leadFieldHeartRegion",instance,mesh,
            IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar("zero",dimless,0));
        heartCells=mesh.cellZones()[zone].size();
        for (label c:mesh.cellZones()[zone]) region[c]=1;
        region.write();
    }
    metadata.add("heartCells",returnReduce(heartCells,sumOp<label>()));
    metadata.add("heartRegionRecorded",props.found("cellZone"));
    const IOdictionary schemes(IOobject("fvSchemes",mesh.time().system(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE,false));
    metadata.add("conductorSchemes",static_cast<const dictionary&>(schemes));
    const IOdictionary solution(IOobject("fvSolution",mesh.time().system(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE,false));
    if (const dictionary* solver=solution.subDict("solvers").findDict("leadPotential"))
        metadata.add("conductorSolver",*solver);
    const word bodyName(props.lookupOrDefault<word>("bathConductivityField", "bodyAndOrgansConductivity"));
    IOobject bodyIO(bodyName, instance, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
    const bool hasBody = props.found("cellZone") && bodyIO.typeHeaderOk<volScalarField>(false);
    metadata.add("bodyConductivityRecorded", hasBody);
    if (hasBody)
    {
        const volScalarField body(IOobject(bodyName, instance, mesh, IOobject::MUST_READ,
            IOobject::NO_WRITE),mesh);
        volScalarField saved(IOobject("leadFieldBodyConductivity",instance,mesh,
            IOobject::NO_READ,IOobject::NO_WRITE),body);
        saved.write();
    }
    if (props.found("leads")) metadata.add("leads", props.subDict("leads"));
    metadata.regIOobject::write();
}

leadFieldData::leadFieldData(const electroStateProvider& state, const dictionary& props, bool heartOnly)
:
    names_(props.subDict("leadFields").toc()), fields_(names_.size()),
    Gi_(), Ge_(), signatures_(meshSignatures(state.mesh()))
{
    const fvMesh& support = state.baseMesh();
    if (names_.empty())
    {
        FatalIOErrorInFunction(props) << "leadFields must contain at least one named field."
            << exit(FatalIOError);
    }
    tmp<volTensorField> gi = readPhysical(support, props, true);
    tmp<volTensorField> ge = readPhysical(support, props, false);
    const IOdictionary metadata(IOobject("leadFieldMetadata", support.time().constant(),
        support, IOobject::MUST_READ, IOobject::NO_WRITE));
    if (metadata.get<label>("formatVersion") != 1
     || metadata.get<word>("representation") != "conductivityWeightedVector")
    {
        FatalErrorInFunction << "Unsupported leadFieldMetadata format/representation."
            << exit(FatalError);
    }
    if (!sameDictionary(props,metadata.subDict("inputProperties")))
    {
        FatalErrorInFunction << "leadFieldProperties differ from the preprocessing inputs. "
            << "Regenerate sensitivities after changing electrodes or conductor settings."
            << exit(FatalError);
    }
    if (support.nGeometricD()!=metadata.get<label>("geometricDimensions"))
        FatalErrorInFunction << "Lead-field geometric dimension/empty-patch mismatch." << exit(FatalError);
    const label cells=returnReduce(support.nCells(),sumOp<label>());
    const label expected=metadata.get<label>(heartOnly ? "heartCells" : "supportCells");
    if (cells!=expected)
    {
        FatalErrorInFunction << "Lead-field support size mismatch. Use importMode heartOnly "
            << "only for an explicitly exported myocardial mesh and fields." << exit(FatalError);
    }
    if (!heartOnly)
    {
        const IOdictionary schemes(IOobject("fvSchemes",support.time().system(),support,
            IOobject::MUST_READ,IOobject::NO_WRITE,false));
        if (!sameDictionary(schemes,metadata.subDict("conductorSchemes")))
            FatalErrorInFunction << "Conductor fvSchemes changed; regenerate lead fields." << exit(FatalError);
        if (metadata.found("conductorSolver"))
        {
            const IOdictionary solution(IOobject("fvSolution",support.time().system(),support,
                IOobject::MUST_READ,IOobject::NO_WRITE,false));
            const dictionary* solver=solution.subDict("solvers").findDict("leadPotential");
            if (!solver || !sameDictionary(*solver,metadata.subDict("conductorSolver")))
                FatalErrorInFunction << "Reciprocal solver controls changed; regenerate lead fields." << exit(FatalError);
        }
    }
    const tmp<volTensorField> savedGi = readStatic<volTensorField>(support, "leadFieldGi", props, sigmaDims);
    const tmp<volTensorField> savedGe = readStatic<volTensorField>(support, "leadFieldGe", props, sigmaDims);
    const tmp<volVectorField> savedSignature = readStatic<volVectorField>(support, "leadFieldMeshSignature", props, dimless);
    const dictionary& fieldDict = props.subDict("leadFields");
    const dictionary& savedFields = metadata.subDict("leadFields");
    if (savedFields.size() != fieldDict.size())
    {
        FatalErrorInFunction << "leadFields do not match leadFieldMetadata. Regenerate sensitivities."
            << exit(FatalError);
    }
    if (props.found("leads") != metadata.found("leads"))
    {
        FatalErrorInFunction << "Lead definitions do not match leadFieldMetadata."
            << exit(FatalError);
    }
    if (props.found("leads"))
    {
        if (!sameDictionary(props.subDict("leads"),metadata.subDict("leads")))
        {
            FatalErrorInFunction << "Lead definitions do not match leadFieldMetadata."
                << exit(FatalError);
        }
    }

    Gi_.reset(new volTensorField(IOobject("leadFieldGiActive", state.mesh().time().timeName(),
        state.mesh(), IOobject::NO_READ, IOobject::NO_WRITE), state.mesh(),
        dimensionedTensor("zero", sigmaDims, Zero)));
    Ge_.reset(new volTensorField(IOobject("leadFieldGeActive", state.mesh().time().timeName(),
        state.mesh(), IOobject::NO_READ, IOobject::NO_WRITE), state.mesh(),
        dimensionedTensor("zero", sigmaDims, Zero)));
    bool provenanceOK = true;
    const vectorField supportSignatures(meshSignatures(support));
    forAll(supportSignatures,c)
        provenanceOK = provenanceOK && supportSignatures[c]==savedSignature()[c];
    if (metadata.get<bool>("heartRegionRecorded"))
    {
        const tmp<volScalarField> region=readStatic<volScalarField>(support,"leadFieldHeartRegion",props,dimless);
        if (heartOnly)
        {
            for (scalar m:region().primitiveField()) provenanceOK=provenanceOK && m==1;
        }
        else
        {
            const label zone=support.cellZones().findZoneID(props.get<word>("cellZone"));
            if (zone<0) FatalErrorInFunction << "Unknown lead-field heart cellZone." << exit(FatalError);
            scalarField expectedMask(support.nCells(),0);
            for (label c:support.cellZones()[zone]) expectedMask[c]=1;
            forAll(expectedMask,c) provenanceOK=provenanceOK && expectedMask[c]==region()[c];
        }
    }
    autoPtr<volScalarField> body, savedBody;
    if (metadata.get<bool>("bodyConductivityRecorded"))
    {
        body.reset(readStatic<volScalarField>(support,
            props.lookupOrDefault<word>("bathConductivityField", "bodyAndOrgansConductivity"),
            props, sigmaDims).ptr());
        savedBody.reset(readStatic<volScalarField>(support,"leadFieldBodyConductivity",
            props,sigmaDims).ptr());
        forAll(body(),c)
            provenanceOK = provenanceOK && std::isfinite(body()[c])
                && mag(body()[c]-savedBody()[c]) <= 1e-12*max(mag(body()[c]),VSMALL);
    }
    forAll(signatures_, c)
    {
        const label s = sourceCell(state, c);
        Gi_()[c] = gi()[s]; Ge_()[c] = ge()[s];
        const scalar scale = max(max(mag(gi()[s]), mag(ge()[s])), VSMALL);
        provenanceOK = provenanceOK && signatures_[c] == savedSignature()[s]
            && mag(gi()[s]-savedGi()[s]) <= 1e-12*scale
            && mag(ge()[s]-savedGe()[s]) <= 1e-12*scale;
    }
    reduce(provenanceOK, andOp<bool>());
    if (!provenanceOK)
    {
        FatalErrorInFunction << "Lead-field mesh or physical conductivity provenance mismatch. "
            << "Regenerate sensitivities for this mesh and Gi/Ge." << exit(FatalError);
    }
    bool finite = true;
    forAll(names_, i)
    {
        const word fieldName(fieldDict.get<word>(names_[i]));
        if (!savedFields.found(names_[i]) || savedFields.get<word>(names_[i]) != fieldName)
        {
            FatalErrorInFunction << "leadFields do not match leadFieldMetadata."
                << exit(FatalError);
        }
        tmp<volVectorField> input = readStatic<volVectorField>(support, fieldName, props, sensitivityDims);
        fields_.set(i, new volVectorField(IOobject("leadFieldActive_"+names_[i],
            state.mesh().time().timeName(), state.mesh(), IOobject::NO_READ,
            IOobject::NO_WRITE), state.mesh(), dimensionedVector("zero", sensitivityDims, Zero)));
        forAll(fields_[i], c)
        {
            fields_[i][c] = input()[sourceCell(state,c)];
            for (direction d=0; d<3; ++d) finite = finite && std::isfinite(fields_[i][c][d]);
        }
    }
    reduce(finite, andOp<bool>());
    if (!finite)
    {
        FatalErrorInFunction << "Lead-field sensitivities contain non-finite values."
            << exit(FatalError);
    }
}

void leadFieldData::project(const electroStateProvider& state, scalarField& values) const
{
    const fvMesh& mesh = state.mesh();
    if (mesh.moving() || mesh.topoChanging() || mesh.nCells()!=signatures_.size())
    {
        FatalErrorInFunction << "Static leadFieldECG requires a fixed mesh."
            << exit(FatalError);
    }
    if (!state.VmPtr())
    {
        FatalErrorInFunction << "leadFieldECG requires Vm." << exit(FatalError);
    }
    if (state.VmPtr()->dimensions() != dimensionSet(1,2,-3,0,0,-1,0))
    {
        FatalErrorInFunction << "leadFieldECG requires Vm in volts." << exit(FatalError);
    }
    // Observation owns its field and wall-state key. Never correct or reuse
    // the source's registered gradient: it can belong to an earlier iterate.
    const volScalarField& source=*state.VmPtr();
    volScalarField voltage(IOobject("leadFieldObservationVm",mesh.time().timeName(),mesh,
        IOobject::NO_READ,IOobject::NO_WRITE),source);
    const labelIOList* sourceStates=mesh.findObject<labelIOList>
        (sealedWallPolicy::statesName(source.name()));
    if (sourceStates)
        sealedWallPolicy::publish(mesh,voltage.name(),
            sealedWallPolicy::states(mesh,source.name()));
    else
        FatalErrorInFunction << "Source " << source.name()
            << " has no published wall-gradient policy. Its owner must publish it."
            << exit(FatalError);
    voltage.correctBoundaryConditions();
    // Select the observation scheme by its configuration entry, but retain
    // the actual field name for wall-policy lookup and gradient cache identity.
    const tmp<volVectorField> observationGradient=fv::gradScheme<scalar>::New
    (
        mesh,mesh.gradScheme("leadFieldGradVm")
    )().grad(voltage,"grad("+voltage.name()+')');
    const vectorField& grad=observationGradient().primitiveField();
    values.setSize(names_.size(), 0.0);
    forAll(names_, i)
    {
        scalar sum = 0;
        forAll(mesh.V(), c) sum += mesh.V()[c]*(fields_[i][c] & grad[c]);
        reduce(sum, sumOp<scalar>());
        values[i] = sum;
    }
}

}
