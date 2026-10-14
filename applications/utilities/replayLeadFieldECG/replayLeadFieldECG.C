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

Application
    replayLeadFieldECG

Description
    Observe stored support-mesh Vm snapshots using reciprocal projection and
    converged passive extracellular recovery. Never advance the cardiac source.
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
#include <cmath>

using namespace Foam;

int main(int argc, char* argv[])
{
    argList::addNote("Compare reciprocal and passive torso ECG on stored Vm.");
    argList::noParallel();
    timeSelector::addOptions();
    argList::addOption("absoluteTolerance", "scalar", "Override passive fixed-point absolute tolerance (V).");
    argList::addOption("relativeTolerance", "scalar", "Override passive fixed-point relative tolerance.");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const instantList times(timeSelector::select0(runTime, args));
    IOdictionary controls(IOobject("passiveECGProperties", runTime.constant(),
        mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    IOdictionary leads(IOobject("leadFieldProperties", runTime.constant(),
        mesh, IOobject::MUST_READ, IOobject::NO_WRITE));
    const dictionary& bathCoeffs = controls.subDict("bathPotentialDomain");
    dictionary heartCoeffs(controls.subDict("bidomainSolverCoeffs"));
    heartCoeffs.add("bathPotentialDomain", bathCoeffs, true);
    if (heartCoeffs.get<word>("cellZone") != leads.get<word>("cellZone")
     || bathCoeffs.getOrDefault<word>("bathConductivityField", "bodyAndOrgansConductivity")
        != leads.getOrDefault<word>("bathConductivityField", "bodyAndOrgansConductivity"))
    {
        FatalIOErrorInFunction(controls)
            << "Forward and reciprocal recording must use the same heart zone "
            << "and passive conductivity field." << exit(FatalIOError);
    }
    if (bathCoeffs.getOrDefault<word>("interfaceConductivityInterpolation",
        "distanceWeightedHarmonic") != "distanceWeightedHarmonic")
    {
        FatalIOErrorInFunction(controls)
            << "Replay currently matches buildLeadFields' distanceWeightedHarmonic "
            << "interface construction. Other conductor schemes need separate controls."
            << exit(FatalIOError);
    }
    for (const word& key : wordList({"groundPatches", "surfaceCurrentPatches"}))
    {
        if (bathCoeffs.found(key) && !bathCoeffs.subDict(key).empty())
        {
            FatalIOErrorInFunction(controls)
                << "Replay matches the insulated passive conductor; " << key
                << " must be empty." << exit(FatalIOError);
        }
    }
    const label maxCorrectors = controls.getOrDefault<label>("maxCorrectors", 200);
    const scalar absoluteTolerance = args.getOrDefault<scalar>("absoluteTolerance",
        controls.getOrDefault<scalar>("absoluteTolerance", 1e-10));
    const scalar relativeTolerance = args.getOrDefault<scalar>("relativeTolerance",
        controls.getOrDefault<scalar>("relativeTolerance", 1e-8));
    const scalar maxLinearResidual = controls.getOrDefault<scalar>("maxLinearResidual", 1e-10);
    if (maxCorrectors < 2 || absoluteTolerance <= 0 || relativeTolerance <= 0
     || maxLinearResidual <= 0 || !std::isfinite(absoluteTolerance)
     || !std::isfinite(relativeTolerance) || !std::isfinite(maxLinearResidual))
    {
        FatalIOErrorInFunction(controls) << "Invalid recovery convergence controls."
            << exit(FatalIOError);
    }

    PtrList<volScalarField> outputFields, postFields;
    autoPtr<ionicModel> ionic;
    autoPtr<electroVerificationModel> verifier;
    autoPtr<myocardiumDomainInterface> heart = myocardiumDomainInterface::New
    (
        mesh, heartCoeffs, outputFields, wordList(), postFields,
        ionic, verifier, runTime.deltaTValue()
    );
    if (verifier.valid())
    {
        FatalIOErrorInFunction(controls)
            << "Passive replay excludes manufactured source forcing."
            << exit(FatalIOError);
    }
    electroVolumeFieldDomain& volume = dynamic_cast<electroVolumeFieldDomain&>(heart());
    volScalarField& Vm = volume.VmRef();
    const labelUList* mapping = heart->subsetCellMapPtr();
    if (!mapping || mapping->size() != Vm.size())
    {
        FatalErrorInFunction << "Replay requires an embedded myocardium cellZone."
            << exit(FatalError);
    }
    extracellularPotentialDomain recovery
    (
        mesh, heart(), heartCoeffs.subDict("bathPotentialDomain")
    );
    forAll(mesh.boundary(), patchI)
    {
        if (!polyPatch::constraintType(mesh.boundary()[patchI].type())
         && recovery.phiE().boundaryField()[patchI].type() != "zeroGradient")
        {
            FatalErrorInFunction << "Passive exterior patch "
                << mesh.boundary()[patchI].name()
                << " must use zeroGradient, matching buildLeadFields."
                << exit(FatalError);
        }
    }
    leadFieldData projection(heart(), leads);
    dictionary torsoDict;
    torsoDict.add("ecgSolver", word("torsoECG"));
    torsoDict.add("electrodePositions", leads.subDict("electrodePositions"));
    ecgDomain torso(recovery, torsoDict, "passiveReplayTorso");
    torsoECG sampler(torsoDict);

    const wordList& names = projection.names();
    scalarField errorSquared(names.size(), 0), referenceSquared(names.size(), 0);
    scalarField maximumError(names.size(), 0);
    const fileName outputDir(runTime.globalPath()/"postProcessing");
    mkDir(outputDir);
    OFstream output(outputDir/"matchedLeadFieldECG.dat");
    OFstream electrodeOutput(outputDir/"matchedTorsoPotentials.dat");
    output.precision(17);
    electrodeOutput.precision(17);
    electrodeOutput << "# time";
    for (const word& name : torso.electrodeNames()) electrodeOutput << ' ' << name;
    electrodeOutput << nl;
    Info << "Passive convergence controls: maxCorrectors=" << maxCorrectors
        << ", absoluteTolerance=" << absoluteTolerance
        << ", relativeTolerance=" << relativeTolerance
        << ", maxLinearResidual=" << maxLinearResidual << nl;
    output << "# time";
    for (const word& name : names)
    {
        output << ' ' << name << "_leadField " << name << "_torso " << name << "_difference";
    }
    output << " correctors maxPhiChange linearResidual sourceUnchanged" << nl;

    label samples = 0;
    forAll(times, timeI)
    {
        runTime.setTime(times[timeI], timeI);
        IOobject sourceIO("Vm", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE, false);
        if (!isFile(sourceIO.objectPath()) && !isFile(sourceIO.objectPath()+".gz"))
        {
            Info << "Skipping " << runTime.timeName() << ": no stored Vm." << nl;
            continue;
        }
        if (!sourceIO.typeHeaderOk<volScalarField>(true))
        {
            FatalErrorInFunction << "Stored Vm has an invalid field header at "
                << runTime.timeName() << exit(FatalError);
        }
        const volScalarField source(sourceIO, mesh);
        if (source.dimensions() != dimensionSet(1, 2, -3, 0, 0, -1, 0))
        {
            FatalErrorInFunction << "Stored Vm must be in volts." << exit(FatalError);
        }
        forAll(Vm, cellI)
        {
            const scalar value = source[(*mapping)[cellI]];
            if (!std::isfinite(value))
            {
                FatalErrorInFunction << "Non-finite stored Vm at "
                    << runTime.timeName() << ", myocardial cell " << cellI
                    << exit(FatalError);
            }
            Vm[cellI] = value;
        }
        Vm.correctBoundaryConditions();
        const scalarField original(Vm.primitiveField());
        List<scalarField> originalBoundary(Vm.boundaryField().size());
        forAll(originalBoundary, patchI) originalBoundary[patchI] = Vm.boundaryField()[patchI];

        recovery.prepareTimeStep(runTime.value(), 0);
        bool converged = false;
        scalar change = GREAT;
        scalar linearResidual = 0;
        label corrector = 0;
        for (; corrector < maxCorrectors; ++corrector)
        {
            const scalarField previous(recovery.phiE().primitiveField());
            recovery.solvePhiE();
            const dictionary& performanceDict = mesh.data().solverPerformanceDict();
            const List<SolverPerformance<scalar>> performances(performanceDict.lookup("phiE"));
            if (performances.empty()
             || !std::isfinite(performances.last().finalResidual())
             || performances.last().finalResidual() > maxLinearResidual)
            {
                FatalErrorInFunction << "Passive phiE linear residual exceeds "
                    << maxLinearResidual << " at " << runTime.timeName()
                    << exit(FatalError);
            }
            linearResidual = max(linearResidual, performances.last().finalResidual());
            change = 0;
            scalar scale = 0;
            forAll(previous, cellI)
            {
                change = max(change, mag(previous[cellI] - recovery.phiE()[cellI]));
                scale = max(scale, mag(recovery.phiE()[cellI]));
            }
            if (corrector > 0 && change < absoluteTolerance + relativeTolerance*scale)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
        {
            FatalErrorInFunction << "Passive recovery did not converge at "
                << runTime.timeName() << "; last max change " << change
                << ". Preserve this failure; do not interpret unconverged traces."
                << exit(FatalError);
        }
        scalarField sampled, projected;
        sampler.solve(torso, runTime.value(), 0, sampled);
        projection.project(heart(), projected);
        if (original != Vm.primitiveField())
        {
            FatalErrorInFunction << "Recording changed frozen Vm." << exit(FatalError);
        }
        forAll(originalBoundary, patchI)
        {
            if (originalBoundary[patchI] != Vm.boundaryField()[patchI])
            {
                FatalErrorInFunction << "Recording changed frozen Vm boundary values."
                    << exit(FatalError);
            }
        }
        electrodeOutput << runTime.value();
        for (const scalar value : sampled) electrodeOutput << ' ' << value;
        electrodeOutput << endl;
        output << runTime.value();
        forAll(names, leadI)
        {
            const dictionary& weights = leads.subDict("leads").subDict(names[leadI]);
            scalar observed = 0;
            forAll(torso.electrodeNames(), electrodeI)
            {
                observed += weights.getOrDefault<scalar>(torso.electrodeNames()[electrodeI], 0)
                    *sampled[electrodeI];
            }
            const scalar difference = projected[leadI] - observed;
            errorSquared[leadI] += sqr(difference);
            referenceSquared[leadI] += sqr(observed);
            maximumError[leadI] = max(maximumError[leadI], mag(difference));
            output << ' ' << projected[leadI] << ' ' << observed << ' ' << difference;
        }
        output << ' ' << corrector + 1 << ' ' << change << ' ' << linearResidual << " 1" << endl;
        ++samples;
        Info << "Passive replay t=" << runTime.timeName() << " converged in "
            << corrector + 1 << " corrections; source unchanged." << nl;
    }
    if (!samples)
    {
        FatalErrorInFunction << "No stored Vm samples selected." << exit(FatalError);
    }
    OFstream summary(outputDir/"matchedLeadFieldECGSummary.dat");
    summary.precision(17);
    summary << "# lead RMS_difference_V RMS_forward_V maximum_difference_V" << nl;
    forAll(names, leadI)
    {
        summary << names[leadI] << ' ' << Foam::sqrt(errorSquared[leadI]/samples)
            << ' ' << Foam::sqrt(referenceSquared[leadI]/samples)
            << ' ' << maximumError[leadI] << nl;
    }
    summary << "samples " << samples << nl << "sourceUnchanged 1" << nl;
    Info << "Completed " << samples << " matched-source samples." << nl << "End" << nl;
    return 0;
}
