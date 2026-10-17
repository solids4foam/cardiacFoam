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
#include "myocardiumDomainInterface.H"
#include "electroVolumeFieldDomain.H"
#include "extracellularPotentialDomain.H"
#include "ecgDomain.H"
#include "torsoECG.H"
#include "leadFieldData.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
using namespace Foam;

int main(int argc,char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    IOstream::defaultPrecision(17);
    IOdictionary electro(IOobject("electroProperties",runTime.constant(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE));
    if (electro.get<word>("myocardiumSolver")!="bidomainSolver")
        FatalErrorInFunction << "Recovery fixture requires bidomainSolver physical Gi/Ge."
            << exit(FatalError);
    IOdictionary props(IOobject("leadFieldProperties",runTime.constant(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE));
    PtrList<volScalarField> outFields,postFields;
    autoPtr<ionicModel> ionic;
    autoPtr<electroVerificationModel> verifier;
    autoPtr<myocardiumDomainInterface> heart=myocardiumDomainInterface::New
        (mesh,electro.subDict("bidomainSolverCoeffs"),outFields,wordList(),
            postFields,ionic,verifier,runTime.deltaTValue());
    if (verifier.valid())
        FatalErrorInFunction << "Passive comparison excludes manufactured source forcing."
            << exit(FatalError);
    electroVolumeFieldDomain& volume=dynamic_cast<electroVolumeFieldDomain&>(heart());
    volScalarField& Vm=volume.VmRef();
    // Prescribed voltage excites an observable nonzero signal. Its scale is
    // artificial; no ionic step or time advance is performed.
    forAll(Vm,c) Vm[c]=0.01*Vm.mesh().C()[c].x();
    Vm.correctBoundaryConditions();
    const scalarField original(Vm.primitiveField());
    List<scalarField> boundaries(Vm.boundaryField().size());
    forAll(boundaries,p) boundaries[p]=Vm.boundaryField()[p];
    dictionary bath;
    bath.add("bathCellZones",wordList(1,"bath"));
    bath.add("bathConductivityField",props.get<word>("bathConductivityField"));
    bath.add("interfaceConductivityInterpolation",word("distanceWeightedHarmonic"));
    bath.add("phiERefPoint",props.get<point>("phiERefPoint"));
    bath.add("phiEReferenceValue",scalar(0));
    dictionary recoveryCoeffs(electro.subDict("bidomainSolverCoeffs"));
    recoveryCoeffs.add("bathPotentialDomain",bath);
    extracellularPotentialDomain recovery(mesh,heart(),
        recoveryCoeffs.subDict("bathPotentialDomain"));
    recovery.prepareTimeStep(runTime.value(),0);
    // Fixed-point iterations of the production corrected assembly with the
    // same frozen Vm. This is algebraic convergence, not time integration.
    bool converged=false;
    for(label i=0;i<100;++i)
    {
        const scalarField old(recovery.phiE().primitiveField());
        recovery.solvePhiE();
        scalar change=0,scale=0;
        forAll(old,c)
        {
            change=max(change,mag(old[c]-recovery.phiE()[c]));
            scale=max(scale,mag(recovery.phiE()[c]));
        }
        reduce(change,maxOp<scalar>()); reduce(scale,maxOp<scalar>());
        if(i>0 && change<1e-11*max(scale,VSMALL)) {converged=true;break;}
    }
    if(!converged) FatalErrorInFunction << "Frozen recovery did not converge." << exit(FatalError);
    dictionary torsoDict;
    torsoDict.add("ecgSolver",word("torsoECG"));
    torsoDict.add("electrodePositions",props.subDict("electrodePositions"));
    ecgDomain torso(recovery,torsoDict,"recoveryTorso");
    torsoECG kernel(torsoDict);
    scalarField sampled;
    kernel.solve(torso,runTime.value(),0,sampled);
    leadFieldData lead(heart(),props);
    scalarField projected;
    lead.project(heart(),projected);
    const wordList channels(props.subDict("leadFields").toc());
    scalar worst=0;
    scalarField recovered(channels.size()),errors(channels.size());
    // Independent discrete oracle for the committed 20-cell heart column:
    // Gi/(Gi+Ge)=1/2 and Vm endpoint-cell values .00025 and .00975 V.
    const scalar expected=0.00475;
    forAll(channels,l)
    {
        scalar observed=0;
        const dictionary& weights=props.subDict("leads").subDict(channels[l]);
        forAll(torso.electrodeNames(),e)
            if(weights.found(torso.electrodeNames()[e]))
                observed+=weights.get<scalar>(torso.electrodeNames()[e])*sampled[e];
        recovered[l]=observed;
        const scalar error=mag(projected[l]-observed);
        errors[l]=error;
        worst=max(worst,error);
        Info<< "channel " << channels[l] << " leadField " << projected[l]
            << " torso " << observed << " absoluteError " << error << nl;
        if(mag(observed-expected)>1e-8 || mag(projected[l]-expected)>1e-8 || error>1e-8)
            FatalErrorInFunction << "Matched passive column mismatch." << exit(FatalError);
    }
    if(original!=Vm.primitiveField())
        FatalErrorInFunction << "Passive recovery changed prescribed Vm." << exit(FatalError);
    forAll(boundaries,p) if(boundaries[p]!=Vm.boundaryField()[p])
        FatalErrorInFunction << "Passive recovery changed source boundaries." << exit(FatalError);
    if(Pstream::master())
    {
        OFstream output(runTime.globalPath()/"postProcessing"/"leadFieldRecoveryCheck.dat");
        output.precision(17);
        output<< "# channel projected torso expected absoluteDifference" << nl;
        forAll(channels,l) output<< channels[l] << ' ' << projected[l] << ' '
            << recovered[l] << ' ' << expected << ' ' << errors[l] << nl;
        output<< "maxAbsoluteError " << worst << nl << "sourceUnchanged 1" << nl;
    }
    Info<< "PASS: frozen-source leadFieldECG/torsoECG comparison" << nl;
    Info<< "End" << nl;
    return 0;
}
