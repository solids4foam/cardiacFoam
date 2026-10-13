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

#include "fvCFD.H"
#include "leadFieldData.H"
#include "extracellularFaceConductivity.H"
#include "PstreamReduceOps.H"
#include "sealedWallPolicy.H"
#include "SVD.H"

#include <cmath>
using namespace Foam;

// A cell observation and a cell reciprocal current use exactly the same
// location/ownership. Nearest-cell approximation is opt-in and distance-limited.
label electrodeCell(const fvMesh& mesh, const point& p, scalar maxDistance)
{
    // Cell-centre construction can communicate on coupled patches. Build
    // it collectively before processor-specific electrode ownership branches.
    const vectorField& centres=mesh.C().primitiveField();
    label cell=mesh.findCell(p);
    label owners=cell>=0 ? 1 : 0;
    reduce(owners,sumOp<label>());
    scalar dist=GREAT;
    if (cell>=0) dist=magSqr(centres[cell]-p);
    if (!owners && maxDistance>0)
    {
        forAll(centres,i)
        {
            const scalar d=magSqr(centres[i]-p);
            if (d<dist) {dist=d;cell=i;}
        }
    }
    scalar nearest=dist;
    reduce(nearest,minOp<scalar>());
    label proc=(cell>=0 && dist==nearest) ? Pstream::myProcNo() : Pstream::nProcs();
    reduce(proc,minOp<label>());
    if ((!owners && (maxDistance<=0 || nearest>sqr(maxDistance))) || proc==Pstream::nProcs())
        FatalErrorInFunction << "Electrode " << p << " has no cell ownership. "
            << "Use an in-domain location or explicitly set maximumElectrodeDistance."
            << exit(FatalError);
    return Pstream::myProcNo()==proc ? cell : -1;
}

// The normal matrix is PSD. Do not build a pseudoinverse from OpenFOAM's
// cubic eigenvalues/eigenvectors: repeated rotated eigenpairs can produce
// non-orthogonal directions and false rank. SVD resolves rank independently;
// use the ordinary inverse in full rank and its range pseudoinverse otherwise.
vector reconstructGradient(const symmTensor& normal,const vector& rhs,
    const label requiredRank,const label cell=-1)
{
    if (tr(normal)<=VSMALL)
        FatalErrorInFunction << "Cannot reconstruct a lead field in an isolated myocardial cell."
            << exit(FatalError);
    scalarRectangularMatrix matrix(3,3);
    matrix(0,0)=normal.xx(); matrix(0,1)=normal.xy(); matrix(0,2)=normal.xz();
    matrix(1,0)=normal.xy(); matrix(1,1)=normal.yy(); matrix(1,2)=normal.yz();
    matrix(2,0)=normal.xz(); matrix(2,1)=normal.yz(); matrix(2,2)=normal.zz();
    const SVD decomposition(matrix,1e-12);
    if(!decomposition.converged())
        FatalErrorInFunction << "Myocardial stencil SVD did not converge at cell "
            << cell << exit(FatalError);
    const label rank=3-decomposition.nZeros();
    if(rank<requiredRank)
        FatalErrorInFunction << "Insufficient myocardial neighbour stencil at cell "
            << cell << ": cannot reconstruct all active gradient directions. "
            << "Refine the heart mesh or import externally reconstructed sensitivities."
            << exit(FatalError);
    if(rank==3) return inv(normal)&rhs;
    const scalarRectangularMatrix inverse(decomposition.VSinvUt());
    vector result(Zero);
    for(direction i=0;i<3;++i) for(direction j=0;j<3;++j)
        result[i]+=inverse(i,j)*rhs[j];
    return result;
}

void checkReconstruction()
{
    const vector axis(vector(1,2,3)/Foam::sqrt(14.0));
    const vector truth(3,-1,2);
    const symmTensor matrices[]={symmTensor::I,2*symmTensor::I+sqr(axis),
        symmTensor::I-sqr(axis),sqr(axis)};
    const label ranks[]={3,3,2,1};
    const vector expected[]={truth,truth,truth-axis*(axis&truth),axis*(axis&truth)};
    for(label i=0;i<4;++i)
    {
        const vector actual=reconstructGradient(matrices[i],matrices[i]&truth,ranks[i]);
        if(mag(actual-expected[i])>1e-10)
            FatalErrorInFunction << "Myocardial affine reconstruction mismatch for rank "
                << ranks[i] << ": " << actual << " versus " << expected[i]
                << exit(FatalError);
        Info<< "PASS: myocardial reconstruction algebra control " << i
            << ", rank " << ranks[i] << ", error " << mag(actual-expected[i]) << nl;
    }
}

// Reconstruct the myocardial side of z using only myocardial cell neighbours.
// Do not differentiate through a conductivity jump, or impose insulation on
// the artificial restriction boundary: reciprocal current crosses it.
vectorField myocardialGradient(const fvMesh& mesh, const volScalarField& z,
    const volScalarField& mask)
{
    labelList active(mesh.nCells(),-1);
    label count=0;
    forAll(active,c) if (mask[c]>0.5) active[c]=count++;
    Field<symmTensor> normal(count,symmTensor::zero);
    vectorField rhs(count,vector::zero);
    auto add=[&](label c, const vector& d, scalar difference)
    {
        const scalar length=magSqr(d);
        if (length>VSMALL)
        {
            normal[active[c]]+=sqr(d)/length;
            rhs[active[c]]+=d*(difference/length);
        }
    };
    const labelUList& own=mesh.owner();
    const labelUList& nei=mesh.neighbour();
    forAll(own,f) if (active[own[f]]>=0 && active[nei[f]]>=0)
    {
        const vector d=mesh.C()[nei[f]]-mesh.C()[own[f]];
        const scalar dz=z[nei[f]]-z[own[f]];
        add(own[f],d,dz); add(nei[f],-d,-dz);
    }
    forAll(mesh.boundary(),p) if (mesh.boundary()[p].coupled())
    {
        const scalarField neighborZ(z.boundaryField()[p].patchNeighbourField());
        const scalarField neighborMask(mask.boundaryField()[p].patchNeighbourField());
        const vectorField delta(mesh.boundary()[p].delta());
        const labelUList& faceCells=mesh.boundary()[p].faceCells();
        forAll(faceCells,f) if (active[faceCells[f]]>=0 && neighborMask[f]>0.5)
            add(faceCells[f],delta[f],neighborZ[f]-z[faceCells[f]]);
    }
    vectorField result(mesh.nCells(),vector::zero);
    forAll(active,c) if (active[c]>=0)
    {
        result[c]=reconstructGradient(normal[active[c]],rhs[active[c]],
            mesh.nGeometricD(),c);
    }
    return result;
}

int main(int argc, char* argv[])
{
    argList::addBoolOption("checkReconstruction","Check repeated-eigenvalue and reduced-rank affine reconstruction algebra.");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    IOstream::defaultPrecision(17);
    if(args.found("checkReconstruction"))
    {
        checkReconstruction();
        Info<< "End" << nl;
        return 0;
    }
    const IOdictionary props(IOobject("leadFieldProperties",runTime.constant(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE));
    const word instance(props.lookupOrDefault<word>("fieldInstance","0"));
    tmp<volTensorField> Gi=leadFieldData::readPhysical(mesh,props,true);
    tmp<volTensorField> Ge=leadFieldData::readPhysical(mesh,props,false);
    const dimensionSet sigmaDims(-1,-3,3,0,0,2,0);
    volScalarField mask(IOobject("leadFieldHeartMask",instance,mesh,IOobject::NO_READ,
        IOobject::NO_WRITE),mesh,dimensionedScalar("one",dimless,1));
    if (props.found("cellZone"))
    {
        mask=dimensionedScalar("zero",dimless,0);
        const word zoneName(props.get<word>("cellZone"));
        const label zone=mesh.cellZones().findZoneID(zoneName);
        if (zone<0) FatalErrorInFunction << "Unknown heart cellZone " << zoneName << exit(FatalError);
        for (label c:mesh.cellZones()[zone]) mask[c]=1;
    }
    label heartCells=0;
    for (scalar m:mask.primitiveField()) if (m>0.5) ++heartCells;
    reduce(heartCells,sumOp<label>());
    if (!heartCells) FatalErrorInFunction << "The lead-field heart region is empty." << exit(FatalError);
    const bool mixed=heartCells<returnReduce(mesh.nCells(),sumOp<label>());
    if (!mixed)
        FatalErrorInFunction << "buildLeadFields requires myocardium embedded in a passive conductor."
            << exit(FatalError);
    forAll(mesh.boundary(),p)
    {
        if (!polyPatch::constraintType(mesh.boundary()[p].type()))
            for (label c:mesh.boundary()[p].faceCells()) if (mask[c]>0.5)
                FatalErrorInFunction << "Myocardium touches physical exterior patch "
                    << mesh.boundary()[p].name() << ". Add a surrounding conductor." << exit(FatalError);
    }
    volTensorField giGlobal(IOobject("leadFieldGiGlobal",instance,mesh,IOobject::NO_READ,
        IOobject::NO_WRITE),mesh,dimensionedTensor("zero",sigmaDims,Zero));
    volTensorField geGlobal(IOobject("leadFieldGeGlobal",instance,mesh,IOobject::NO_READ,
        IOobject::NO_WRITE),Ge());
    if (mixed)
    {
        volScalarField body(IOobject(props.lookupOrDefault<word>("bathConductivityField",
            "bodyAndOrgansConductivity"),instance,mesh,IOobject::MUST_READ,
            IOobject::NO_WRITE),mesh);
        if (body.dimensions()!=sigmaDims)
            FatalErrorInFunction << "Body conductivity must be in S/m." << exit(FatalError);
        forAll(geGlobal,c) if (mask[c]<0.5)
        {
            if (!std::isfinite(body[c]) || body[c]<=0)
                FatalErrorInFunction << "Body conductivity must be finite and positive." << exit(FatalError);
            geGlobal[c]=body[c]*tensor::I;
        }
    }
    forAll(giGlobal,c) if (mask[c]>0.5) giGlobal[c]=Gi()[c];
    mask.correctBoundaryConditions();
    giGlobal.correctBoundaryConditions(); geGlobal.correctBoundaryConditions();
    surfaceTensorField coeff(fvc::interpolate(geGlobal)+fvc::interpolate(giGlobal));
    const labelUList& own=mesh.owner();
    const labelUList& nei=mesh.neighbour();
    const scalarField& weights=mesh.weights().primitiveField();
    forAll(own,f)
    {
        if ((mask[own[f]]>0.5)!=(mask[nei[f]]>0.5))
            coeff[f]=extracellularFaceConductivity::distanceWeightedHarmonic(
                geGlobal[own[f]],geGlobal[nei[f]],1-weights[f],weights[f]);
    }
    forAll(coeff.boundaryField(),p)
    {
        if (mesh.boundary()[p].coupled())
        {
            const scalarField mOwn(mask.boundaryField()[p].patchInternalField());
            const scalarField mNei(mask.boundaryField()[p].patchNeighbourField());
            const tensorField gOwn(geGlobal.boundaryField()[p].patchInternalField());
            const tensorField gNei(geGlobal.boundaryField()[p].patchNeighbourField());
            const scalarField& w=mesh.weights().boundaryField()[p];
            forAll(mOwn,f) if ((mOwn[f]>0.5)!=(mNei[f]>0.5))
                coeff.boundaryFieldRef()[p][f]=extracellularFaceConductivity::distanceWeightedHarmonic(
                    gOwn[f],gNei[f],1-w[f],w[f]);
        }
        else coeff.boundaryFieldRef()[p]=geGlobal.boundaryField()[p].patchInternalField();
    }
    wordList types(mesh.boundary().size(),"zeroGradient");
    forAll(types,p) if (mesh.boundary()[p].coupled() || mesh.boundary()[p].type()=="empty")
        types[p]=mesh.boundary()[p].type();
    const dimensionSet resistance(1,2,-3,0,0,-2,0);
    const dictionary& positions=props.subDict("electrodePositions");
    const wordList electrodes(positions.toc());
    labelList cells(electrodes.size());
    forAll(electrodes,i) cells[i]=electrodeCell(mesh,positions.get<point>(electrodes[i]),
        props.lookupOrDefault<scalar>("maximumElectrodeDistance",0));
    const point refPoint(props.get<point>("phiERefPoint"));
    const label ref=electrodeCell(mesh,refPoint,0);
    const dictionary& leads=props.subDict("leads");
    const dictionary& fields=props.subDict("leadFields");
    const wordList names(fields.toc());
    if (names.empty() || names.size()!=leads.size())
        FatalErrorInFunction << "leadFields and leads must have matching non-empty channel sets." << exit(FatalError);
    for (const word& name:names)
    {
        const dictionary& lead=leads.subDict(name);
        scalar balance=0, strength=0;
        volScalarField source(IOobject("leadFieldCurrent",instance,mesh,IOobject::NO_READ,
            IOobject::NO_WRITE),mesh,dimensionedScalar("zero",dimensionSet(0,-3,0,0,0,0,0),0),types);
        forAllConstIter(dictionary,lead,iter)
        {
            const word& e=iter().keyword();
            const scalar weight=lead.get<scalar>(e);
            if (!std::isfinite(weight)) FatalErrorInFunction << "Non-finite lead weight." << exit(FatalError);
            label index=-1;
            forAll(electrodes,i) if (electrodes[i]==e) index=i;
            if (index<0) FatalErrorInFunction << "Unknown electrode " << e << exit(FatalError);
            balance+=weight; strength+=mag(weight);
            if (cells[index]>=0) source[cells[index]]+=weight/mesh.V()[cells[index]];
        }
        if (!strength || mag(balance)>1e-12*strength)
            FatalErrorInFunction << "Lead " << name << " must have nonzero balanced electrode weights." << exit(FatalError);
        volScalarField z(IOobject("leadPotential_"+name,instance,mesh,IOobject::NO_READ,
            IOobject::NO_WRITE),mesh,dimensionedScalar("zero",resistance,0),
            types);
        volVectorField gradZ(IOobject("grad("+z.name()+")",instance,mesh,
            IOobject::NO_READ,IOobject::NO_WRITE),mesh,
            dimensionedVector("zero",resistance/dimLength,Zero));
        sealedWallPolicy::publish(z,true);
        scalar previous=GREAT;
        bool converged=false;
        const label maxCorrectors=props.lookupOrDefault<label>("maxCorrectors",100);
        const scalar correctionTolerance=props.lookupOrDefault<scalar>("correctionTolerance",1e-10);
        if (maxCorrectors<2 || correctionTolerance<=0)
            FatalErrorInFunction << "maxCorrectors must be >=2 and correctionTolerance positive." << exit(FatalError);
        for (label corr=0;corr<maxCorrectors;++corr)
        {
            gradZ=fvc::grad(z);
            z.correctBoundaryConditions();
            const scalarField old(z.primitiveField());
            fvScalarMatrix eqn(-fvm::laplacian(coeff,z)==source);
            if (ref>=0) eqn.setReference(ref,0,true);
            const solverPerformance performance=eqn.solve("leadPotential");
            if (!performance.converged())
                FatalErrorInFunction << "Reciprocal linear solve did not converge for "
                    << name << exit(FatalError);
            scalar delta=0, magnitude=0;
            forAll(z,c) {delta=max(delta,mag(z[c]-old[c]));magnitude=max(magnitude,mag(z[c]));}
            reduce(delta,maxOp<scalar>()); reduce(magnitude,maxOp<scalar>());
            previous=delta/max(magnitude,VSMALL);
            if (previous<correctionTolerance) {converged=true;break;}
        }
        if (!converged) FatalErrorInFunction << "Lead " << name
            << " did not converge corrected assemblies; relative update " << previous << exit(FatalError);
        const vectorField grad=myocardialGradient(mesh,z,mask);
        volVectorField L(IOobject(fields.get<word>(name),instance,mesh,IOobject::NO_READ,
            IOobject::NO_WRITE),mesh,dimensionedVector("zero",dimensionSet(0,-2,0,0,0,0,0),Zero));
        forAll(L,c) L[c]=-(giGlobal[c].T() & grad[c]);
        L.write();
        if (props.lookupOrDefault<Switch>("writePotentials",false))
        {
            geGlobal.write(); gradZ.write(); z.write();
        }
        Info<< "Built lead " << name << " (positive potential minus weighted reference)." << nl;
    }
    leadFieldData::writeProvenance(mesh,Gi(),Ge(),props,"myocardialNeighbourLeastSquaresSVD");
    Info<< "End" << nl;
    return 0;
}
