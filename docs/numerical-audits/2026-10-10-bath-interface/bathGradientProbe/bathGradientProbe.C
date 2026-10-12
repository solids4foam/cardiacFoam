// Temporary diagnostic: exact fields, exact face traces, and geometry controls.
#include "fvCFD.H"
#include "fvMeshSubset.H"
#include "../../../src/verificationModels/bathBidomainVerification/manufacturedTransmissionBathBidomainReference.H"
using namespace Foam;
namespace M = Foam::transmissionMMS;
bool affineWall = false;
vector wallGradient(const M::Params& p)
{ vector n(Foam::cos(p.alpha),Foam::sin(p.alpha),0),t(-n.y(),n.x(),0); symmTensor G=M::Gi(point(0.9,0.4,0),p); return t-n*((n&G&t)/(n&G&n)); }

wordList types(const fvMesh& m)
{
    wordList out(m.boundary().size(), "fixedValue");
    forAll(out,p) if(polyPatch::constraintType(m.boundary()[p].type())) out[p]=m.boundary()[p].type();
    return out;
}
scalar value(const point& x, const M::Params& p, bool vm)
{
    if(vm && affineWall)return wallGradient(p)&x;
    return vm ? M::Vm(x,0,p) : M::phiE(x,0,M::zoneOf(x,p),p);
}
vector exactGrad(const point& x,const M::Params& p,bool vm)
{
    if(vm && affineWall)return wallGradient(p);
    return vm ? M::gradVm(x,0,p) : M::gradPhiE(x,0,M::zoneOf(x,p),p);
}
void fill(volScalarField& f,const M::Params& p,bool vm)
{
    const fvMesh& m=f.mesh();
    forAll(f,c) f[c]=value(m.C()[c],p,vm);
    forAll(m.boundary(),b)
    {
        if(m.boundary()[b].coupled() || isA<emptyPolyPatch>(m.boundaryMesh()[b]))continue;
        forAll(m.boundary()[b],i)
        {
            const point x=m.Cf().boundaryField()[b][i];
            f.boundaryFieldRef()[b][i]=vm ? value(x,p,true) : M::phiE(x,0,M::myocardium,p);
        }
    }
    f.correctBoundaryConditions();
}
vectorField fit(const fvMesh& m,const M::Params& p,bool vm,bool separate,bool projected)
{
    symmTensorField D(m.nCells(),Zero);vectorField R(m.nCells(),Zero);
    const auto add=[&](label c,vector d,scalar du,scalar area)
    { scalar w=area/magSqr(d);D[c]+=w*sqr(d);R[c]+=w*du*d;};
    const auto& C=m.C();
    forAll(m.neighbour(),f)
    {
        label o=m.owner()[f],n=m.neighbour()[f];
        bool jump=separate && M::zoneOf(C[o],p)!=M::zoneOf(C[n],p);
        if(jump)
        {
            point x=m.Cf()[f]; scalar u=M::phiE(x,0,M::myocardium,p);
            add(o,x-C[o],u-value(C[o],p,vm),m.magSf()[f]);
            add(n,x-C[n],u-value(C[n],p,vm),m.magSf()[f]);
        }
        else
        {
            vector d=C[n]-C[o];scalar du=value(C[n],p,vm)-value(C[o],p,vm);
            add(o,d,du,m.magSf()[f]);add(n,-d,-du,m.magSf()[f]);
        }
    }
    forAll(m.boundary(),b)
    {
        const fvPatch& patch=m.boundary()[b];
        if(isA<emptyPolyPatch>(m.boundaryMesh()[b]))continue;
        vectorField nf(patch.nf()),delta(patch.delta());
        forAll(patch,i)
        {
            label c=patch.faceCells()[i];point x=m.Cf().boundaryField()[b][i];
            vector d;
            scalar u;
            if(patch.coupled())
            {d=delta[i];u=value(C[c]+d,p,vm);}
            else
            {d=x-C[c]; if(projected)d=nf[i]*(nf[i]&d);u=vm ? value(x,p,true):M::phiE(x,0,M::zoneOf(C[c],p),p);}
            add(c,d,u-value(C[c],p,vm),m.magSf().boundaryField()[b][i]);
        }
    }
    vectorField g(m.nCells(),Zero);
    forAll(g,c){D[c].zz()=1;g[c]=inv(D[c])&R[c];}
    return g;
}
void report(const word& labelName,const fvMesh& m,const vectorField& g,const M::Params& p,bool vm,const boolList& wall)
{
    scalar e2=0,t2=0,em=0,v=0,ei2=0,vi=0;
    vector n(Foam::cos(p.alpha),Foam::sin(p.alpha),0);
    boolList exclude(m.nCells(),false);
    if(affineWall)
    {
        forAll(m.boundary(),b)if(m.boundary()[b].coupled())
            for(label c:m.boundary()[b].faceCells())exclude[c]=true;
        for(label iter=0;iter<2;++iter)
        {
            boolList next(exclude);
            forAll(m.neighbour(),f){label o=m.owner()[f],nb=m.neighbour()[f];if(exclude[o]||exclude[nb]){next[o]=true;next[nb]=true;}}
            exclude=next;
        }
    }
    forAll(g,c)
    {
        if(exclude[c])continue;
        vector e=g[c]-exactGrad(m.C()[c],p,vm);scalar V=m.V()[c];
        if(wall[c]){e2+=V*magSqr(e);t2+=V*magSqr(e-n*(n&e));em=max(em,mag(e));v+=V;}
        else{ei2+=V*magSqr(e);vi+=V;}
    }
    Info<<"PROBE "<<labelName<<" wallL2="<<Foam::sqrt(e2/max(v,SMALL))<<" wallMax="<<em<<" tangentL2="<<Foam::sqrt(t2/max(v,SMALL))<<" innerL2="<<Foam::sqrt(ei2/max(vi,SMALL))<<nl;
}
int main(int argc,char* argv[])
{
    argList::noParallel();argList::addBoolOption("linearWall", "affine conormal-null heart field");argList::addBoolOption("affine","affine exact solution");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    M::Params p=M::defaultParams(); affineWall=args.found("linearWall"); if(args.found("affine"))p.solution=0;
    boolList wall(mesh.nCells(),false);
    forAll(mesh.neighbour(),f)
    {label o=mesh.owner()[f],n=mesh.neighbour()[f];if(M::zoneOf(mesh.C()[o],p)!=M::zoneOf(mesh.C()[n],p)){wall[o]=true;wall[n]=true;}}
    volScalarField phi(IOobject("probePhi",runTime.timeName(),mesh),mesh,dimensionedScalar("zero",dimless,0),types(mesh));
    fill(phi,p,false);
    // Outer Dirichlet faces are bath faces, unlike heart trace faces.
    forAll(mesh.boundary(),b)if(!mesh.boundary()[b].coupled()&&!isA<emptyPolyPatch>(mesh.boundaryMesh()[b]))
        forAll(mesh.boundary()[b],i)phi.boundaryFieldRef()[b][i]=value(mesh.Cf().boundaryField()[b][i],p,false);
    phi.correctBoundaryConditions();
    report("baseStock",mesh,fvc::grad(phi)().primitiveField(),p,false,wall);
    report("baseExactJumpTrace",mesh,fit(mesh,p,false,true,false),p,false,wall);
    fvMeshSubset sub(mesh);
    label z=mesh.cellZones().findZoneID("myocardium"); label bp=mesh.boundaryMesh().findPatchID("exposed");
    sub.setCellSubset(mesh.cellZones()[z],bp);
    const fvMesh& heart=sub.subMesh();boolList hw(heart.nCells(),false);
    forAll(heart.boundary(),b)if(!heart.boundary()[b].coupled()&&!isA<emptyPolyPatch>(heart.boundaryMesh()[b]))
        for(label c:heart.boundary()[b].faceCells())hw[c]=true;
    volScalarField vm(IOobject("probeVm",runTime.timeName(),heart),heart,dimensionedScalar("zero",dimless,0),types(heart));fill(vm,p,true);
    report("heartStockExactTrace",heart,fvc::grad(vm)().primitiveField(),p,true,hw);
    report("heartFitTrueDelta",heart,fit(heart,p,true,false,false),p,true,hw);
    report("heartFitProjectedDelta",heart,fit(heart,p,true,false,true),p,true,hw);
    if (!affineWall)
    {
        volScalarField pe(IOobject("probeHeartPhiE",runTime.timeName(),heart),heart,dimensionedScalar("zero",dimless,0),types(heart));fill(pe,p,false);
        report("heartPhiEStockExactTrace",heart,fvc::grad(pe)().primitiveField(),p,false,hw);
        report("heartPhiEFitTrueDelta",heart,fit(heart,p,false,false,false),p,false,hw);
        report("heartPhiEFitProjectedDelta",heart,fit(heart,p,false,false,true),p,false,hw);
    }
    return 0;
}
