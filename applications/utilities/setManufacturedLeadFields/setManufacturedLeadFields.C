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
#include "IOdictionary.H"
#include "sealedWallPolicy.H"
using namespace Foam;
class affineState : public electroStateProvider
{
    const volScalarField& voltage_;
    const volTensorField& conductivity_;
    const volScalarField* phi_;
public:
    affineState(const volScalarField& voltage,const volTensorField& conductivity,
        const volScalarField* phi) : voltage_(voltage),conductivity_(conductivity),phi_(phi) {}
    const fvMesh& mesh() const override { return voltage_.mesh(); }
    const volScalarField* VmPtr() const override { return &voltage_; }
    const volScalarField* phiEPtr() const override { return phi_; }
    const volTensorField* intracellularConductivityPtr() const override { return &conductivity_; }
};

int main(int argc, char* argv[])
{
    argList::addBoolOption("checkSealedWall","Verify private sealed-wall gradient and source independence with affine voltage.");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    IOstream::defaultPrecision(17);
    const IOdictionary props(IOobject("leadFieldProperties", runTime.constant(),
        mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    const scalarField scales(props.lookup("manufacturedScales"));
    const wordList names(props.subDict("leadFields").toc());
    if (scales.size()!=names.size())
        FatalErrorInFunction << "manufacturedScales must match leadFields." << exit(FatalError);
    tmp<volTensorField> Gi = leadFieldData::readPhysical(mesh,props,true);
    tmp<volTensorField> Ge = leadFieldData::readPhysical(mesh,props,false);
    forAll(names,i)
    {
        volVectorField field(IOobject(props.subDict("leadFields").get<word>(names[i]),
            props.lookupOrDefault<word>("fieldInstance","0"), mesh,
            IOobject::NO_READ,IOobject::NO_WRITE), mesh,
            dimensionedVector("zero",dimensionSet(0,-2,0,0,0,0,0),Zero));
        forAll(field,c)
        {
            const point& p=mesh.C()[c];
            const scalar x=constant::mathematical::pi*p.x();
            const scalar y=2*constant::mathematical::pi*p.y();
            const scalar z=3*constant::mathematical::pi*p.z();
            field[c]=scales[i]*vector(Foam::sin(x)*Foam::cos(y)*Foam::cos(z),
                2*Foam::cos(x)*Foam::sin(y)*Foam::cos(z),3*Foam::cos(x)*Foam::cos(y)*Foam::sin(z));
        }
        if (args.found("checkSealedWall")) field=dimensionedVector("constant",field.dimensions(),vector(scales[i],0,0));
        field.write();
    }
    leadFieldData::writeProvenance(mesh,Gi(),Ge(),props,"analyticManufactured");
    if (args.found("checkSealedWall"))
    {
        const dimensionSet volts(1,2,-3,0,0,-1,0);
        volTensorField modelGi(IOobject("sealedWallTestGi",runTime.timeName(),mesh,
            IOobject::NO_READ,IOobject::NO_WRITE),mesh,
            dimensionedTensor("Gi",dimensionSet(-1,-3,3,0,0,2,0),tensor(2,0.5,0,0.5,2,0,0,0,1)));
        for (label offset=0;offset<2;++offset)
        {
            wordList types(mesh.boundary().size(),"fixedValue");
            forAll(types,p) if (mesh.boundary()[p].name()=="inlet" || mesh.boundary()[p].name()=="outlet")
                types[p]="zeroGradient";
            volScalarField voltage(IOobject("Vm",runTime.timeName(),mesh,IOobject::NO_READ,
                IOobject::NO_WRITE),mesh,dimensionedScalar("zero",volts,0),types);
            volVectorField stale(IOobject("grad(Vm)",runTime.timeName(),mesh,IOobject::NO_READ,
                IOobject::NO_WRITE),mesh,dimensionedVector("zero",volts/dimLength,Zero));
            volScalarField phi(IOobject("sealedWallTestPhi",runTime.timeName(),mesh,IOobject::NO_READ,
                IOobject::NO_WRITE),mesh,dimensionedScalar("zero",volts,0),
                wordList(mesh.boundary().size(),"fixedValue"));
            const scalar alpha=offset ? 0.4 : 0.25;
            forAll(voltage,c)
            {
                voltage[c]=mesh.C()[c].y()-alpha*mesh.C()[c].x();
                phi[c]=0.1*mesh.C()[c].x()+0.2*mesh.C()[c].y();
            }
            forAll(mesh.boundary(),p)
            {
                scalarField v(mesh.boundary()[p].size()), f(v.size());
                forAll(v,i)
                {
                    const point x=mesh.Cf().boundaryField()[p][i];
                    v[i]=x.y()-alpha*x.x(); f[i]=0.1*x.x()+0.2*x.y();
                }
                voltage.boundaryFieldRef()[p]==v;
                phi.boundaryFieldRef()[p]==f;
            }
            sealedWallPolicy::publish(voltage,true);
            const scalarField original(voltage.primitiveField());
            List<scalarField> boundaries(voltage.boundaryField().size());
            forAll(boundaries,p) boundaries[p]=voltage.boundaryField()[p];
            const scalarField originalPhi(phi.primitiveField());
            List<scalarField> phiBoundaries(phi.boundaryField().size());
            forAll(phiBoundaries,p) phiBoundaries[p]=phi.boundaryField()[p];
            affineState state(voltage,modelGi,offset ? &phi : nullptr);
            leadFieldData data(state,props);
            scalarField values;
            data.project(state,values);
            forAll(values,i) if (mag(values[i]+alpha*scales[i])>1e-9)
                FatalErrorInFunction << "Sealed-wall affine projection mismatch: " << values << exit(FatalError);
            if (voltage.primitiveField()!=original)
                FatalErrorInFunction << "Observation changed source voltage." << exit(FatalError);
            forAll(boundaries,p) if (boundaries[p]!=voltage.boundaryField()[p])
                FatalErrorInFunction << "Observation changed source boundary values." << exit(FatalError);
            if (phi.primitiveField()!=originalPhi)
                FatalErrorInFunction << "Observation changed source extracellular voltage." << exit(FatalError);
            forAll(phiBoundaries,p) if (phiBoundaries[p]!=phi.boundaryField()[p])
                FatalErrorInFunction << "Observation changed source extracellular boundaries." << exit(FatalError);
            Info<< "PASS: private sealed-wall gradient, phiE source snapshot=" << offset << nl;
        }
    }
    Info<< "End" << nl;
    return 0;
}
