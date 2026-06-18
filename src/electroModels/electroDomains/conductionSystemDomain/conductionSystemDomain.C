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

#include "conductionSystemDomain.H"

#include "DynamicList.H"
#include "IOdictionary.H"
#include "PstreamReduceOps.H"
#include "ionicVariableCompatibility.H"

namespace Foam
{

defineTypeNameAndDebug(conductionSystemDomain, 0);

namespace
{

bool selectedSolverRequiresIonicModel(const dictionary& dict)
{
    const word solverType
    (
        dict.lookupOrDefault<word>
        (
            "conductionSystemSolver",
            "monodomain1DSolver"
        )
    );

    return solverType != "eikonalSolver1D"
        && solverType != "eikonalSolver"
        && solverType != "restitutionEikonalSolver1D";
}

} // End anonymous namespace


autoPtr<conductionSystemDomain> conductionSystemDomain::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    scalar initialDeltaT
)
{
    const word domainType
    (
        dict.lookupOrDefault<word>
        (
            "conductionSystemDomain",
            "purkinjeGraphModel"
        )
    );

    if (domainType != "purkinjeGraphModel")
    {
        FatalErrorInFunction
            << "Unsupported conductionSystemDomain type '" << domainType
            << "'. Only 'purkinjeGraphModel' is supported."
            << exit(FatalError);
    }

    return autoPtr<conductionSystemDomain>
    (
        new conductionSystemDomain(mesh, dict, initialDeltaT)
    );
}


void conductionSystemDomain::readGraphFile(const dictionary& dict)
{
    if (!dict.found("graphFile"))
    {
        return;
    }

    const word graphFile(dict.get<word>("graphFile"));

    IOdictionary graphDict
    (
        IOobject
        (
            graphFile,
            supportMesh_.time().caseConstant(),
            supportMesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    graph_.readFromDict(graphDict, reportSetup_);

    const scalar purkinjeConductivity =
        dict.lookupOrDefault<scalar>("purkinjeConductivity", 1.0);

    graph_.edgeConductances *= purkinjeConductivity;

    Info<< "Purkinje edge conductance multiplier: "
        << purkinjeConductivity << nl << endl;

    rootNode_ = graphDict.lookupOrDefault<label>("rootNode", 0);
    terminalNodes_ = labelList(graphDict.lookup("pvjNodes"));
    nodeLocations_ = pointField(graphDict.lookup("points"));
    terminalLocations_ = pointField(graphDict.lookup("pvjLocations"));

    if (rootNode_ < 0 || rootNode_ >= graph_.nNodes)
    {
        FatalErrorInFunction
            << "rootNode " << rootNode_ << " is outside graph node range [0,"
            << graph_.nNodes - 1 << "]."
            << exit(FatalError);
    }

    if (terminalNodes_.size() != terminalLocations_.size())
    {
        FatalErrorInFunction
            << "pvjNodes and pvjLocations from graph file '" << graphFile
            << "' must have the same size. pvjNodes.size()="
            << terminalNodes_.size()
            << " pvjLocations.size()=" << terminalLocations_.size()
            << exit(FatalError);
    }

    if (nodeLocations_.size() != graph_.nNodes)
    {
        FatalErrorInFunction
            << "points from graph file '" << graphFile
            << "' must have one coordinate per graph node. points.size()="
            << nodeLocations_.size()
            << " graph nodes=" << graph_.nNodes
            << exit(FatalError);
    }

    forAll(terminalNodes_, i)
    {
        if (terminalNodes_[i] < 0 || terminalNodes_[i] >= graph_.nNodes)
        {
            FatalErrorInFunction
                << "pvjNodes[" << i << "]=" << terminalNodes_[i]
                << " is outside graph node range [0,"
                << graph_.nNodes - 1 << "]."
                << exit(FatalError);
        }
    }

    if (reportSetup_)
    {
        Info<< "Purkinje graph file '" << graphFile << "': rootNode="
            << rootNode_ << ", terminals=" << terminalNodes_.size() << nl;
        forAll(terminalNodes_, i)
        {
            Info<< "  terminal" << i
                << " node=" << terminalNodes_[i]
                << " location=" << terminalLocations_[i] << nl;
        }
        Info<< endl;
    }
}


void conductionSystemDomain::readRootStimulus(const dictionary& dict)
{
    if (!dict.found("rootStimulus"))
    {
        rootStartTime_ = GREAT;
        rootDuration_ = 0.0;
        rootIntensity_ = 0.0;
        return;
    }

    const dictionary& rsDict = dict.subDict("rootStimulus");

    rootStartTime_ = rsDict.get<scalar>("startTime");
    rootDuration_ = rsDict.lookupOrDefault<scalar>("duration", 0.0);
    rootIntensity_ = rsDict.lookupOrDefault<scalar>("intensity", 0.0);

    if (rsDict.found("node"))
    {
        rootNode_ = rsDict.get<label>("node");
    }

    Info<< "Purkinje root stimulus: node=" << rootNode_
        << ", start=" << rootStartTime_
        << ", duration=" << rootDuration_
        << ", intensity=" << rootIntensity_ << nl << endl;
}


void conductionSystemDomain::initialiseState(const scalar initialDeltaT)
{
    label N = graph_.nNodes;
    label nodesPerProc = N / Pstream::nProcs();
    localStartNode_ = Pstream::myProcNo() * nodesPerProc;
    label endNode = (Pstream::myProcNo() == Pstream::nProcs() - 1) ? N : localStartNode_ + nodesPerProc;
    nLocalNodes_ = endNode - localStartNode_;

    if (selectedSolverRequiresIonicModel(coeffsDict_))
    {
        ionicModelPtr_ = ionicModel::New
        (
            coeffsDict_,
            nLocalNodes_,
            initialDeltaT,
            false
        );
    }

    const scalar vmRest =
    (
        ionicModelPtr_.valid()
      ? ionicModelPtr_->vmRest()
      : coeffsDict_.lookupOrDefault<scalar>("vm1DRest", -0.084)
    );

    Vm1D_.setSize(graph_.nNodes, vmRest);
    Iion1D_.setSize(graph_.nNodes, 0.0);
    activationTime_.setSize(graph_.nNodes, -1.0);

    if (rootStartTime_ < GREAT && rootStartTime_ <= SMALL)
    {
        activationTime_[rootNode_] = rootStartTime_;
    }

    terminalCurrent_.setSize(terminalNodes_.size(), 0.0);
    terminalSource_.setSize(terminalNodes_.size(), 0.0);

    if (ionicModelPtr_.valid())
    {
        Info<< "conductionSystemDomain ionic model: "
            << ionicModelPtr_->type() << nl << endl;
    }
}


void conductionSystemDomain::initialiseOutputControls()
{
    const dictionary& ovDict = coeffsDict_.subOrEmptyDict("outputVariables");

    wordList userExport = ovDict.getOrDefault<wordList>
    (
        "export",
        wordList{"Vm", "IcouplingSource"}
    );

    if (ionicModelPtr_.valid())
    {
        purkinjeModelIO::ResolvedTokens resolved = purkinjeModelIO::filterTokens
        (
            userExport,
            ionicModelPtr_->ioStateNames(),
            ionicModelPtr_->ioNumStates(),
            ionicModelPtr_->ioAlgebraicNames(),
            ionicModelPtr_->ioNumAlgebraic()
        );

        exportVars_ = resolved.networkTokens;
        ionicExport_ = resolved.ionicTokens;

        if (resolved.unknownTokens.size() > 0)
        {
            WarningInFunction
                << "The following export variables are unknown and will be ignored: "
                << resolved.unknownTokens << endl;
        }

        ionicExportStateIndices_.setSize(ionicExport_.size(), -1);
        ionicExportAlgebraicIndices_.setSize(ionicExport_.size(), -1);

        forAll(ionicExport_, i)
        {
            const word& var = ionicExport_[i];
            bool isVmDummy = false;
            label sIdx = -1;
            label aIdx = -1;
            label rIdx = -1;

            ionicVariableCompatibility::resolveVariable
            (
                var,
                ionicModelPtr_->ioStateNames(),
                ionicModelPtr_->ioNumStates(),
                ionicModelPtr_->ioAlgebraicNames(),
                ionicModelPtr_->ioNumAlgebraic(),
                isVmDummy,
                sIdx,
                aIdx,
                rIdx
            );

            if (sIdx >= 0)
            {
                ionicExportStateIndices_[i] = sIdx;
            }
            else if (aIdx >= 0)
            {
                ionicExportAlgebraicIndices_[i] = aIdx;
            }
        }
    }
    else
    {
        exportVars_ = userExport;
    }

    debugVars_ = ovDict.getOrDefault<wordList>
    (
        "debug",
        wordList()
    );

    if (ovDict.found("probeNodes"))
    {
        probeNodes_ = ovDict.lookup("probeNodes");
    }
    else
    {
        probeNodes_ = identity(graph_.nNodes);
    }
}


void conductionSystemDomain::openOutputFile()
{
    if (!Pstream::master())
    {
        return;
    }

    DynamicList<word> colNames;

    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            for (label nodeI = 0; nodeI < graph_.nNodes; ++nodeI)
            {
                colNames.append("node" + Foam::name(nodeI) + "_Vm_V");
            }
        }
        else if (var == "Iion")
        {
            for (label nodeI = 0; nodeI < graph_.nNodes; ++nodeI)
            {
                colNames.append("node" + Foam::name(nodeI) + "_Iion");
            }
        }
        else if (var == "activationTime")
        {
            for (label nodeI = 0; nodeI < graph_.nNodes; ++nodeI)
            {
                colNames.append("node" + Foam::name(nodeI) + "_activationTime");
            }
        }
        else if (var == "IcouplingSource")
        {
            forAll(terminalNodes_, i)
            {
                colNames.append("pvj" + Foam::name(i) + "_IcouplingSource_Am3");
            }
        }
        else if (var == "IcouplingCurrent")
        {
            forAll(terminalNodes_, i)
            {
                colNames.append("pvj" + Foam::name(i) + "_IcouplingCurrent");
            }
        }
    }

    const fileName outDir(time().globalPath()/"postProcessing");
    outputPtr_ = purkinjeModelIO::openTimeSeries
    (
        outDir,
        "purkinjeNetwork.dat",
        colNames
    );

    if (reportSetup_)
    {
        Info<< "conductionSystemDomain: writing to "
            << outDir/"purkinjeNetwork.dat" << nl << endl;
    }

    ionicOutputPtrs_.setSize(ionicExport_.size());
    forAll(ionicExport_, i)
    {
        const word& var = ionicExport_[i];
        DynamicList<word> ionicColNames(probeNodes_.size());
        forAll(probeNodes_, pI)
        {
            ionicColNames.append("node" + Foam::name(probeNodes_[pI]) + "_" + var);
        }

        ionicOutputPtrs_.set
        (
            i,
            purkinjeModelIO::openTimeSeries
            (
                outDir,
                "purkinjeNetwork_" + var + ".dat",
                ionicColNames
            ).ptr()
        );
    }
}


conductionSystemDomain::conductionSystemDomain
(
    const fvMesh& mesh,
    const dictionary& dict,
    const scalar initialDeltaT
)
:
    supportMesh_(mesh),
    coeffsDict_(dict.subDict("purkinjeGraphModelCoeffs")),
    graph_(),
    solverPtr_(conductionSystemSolver::New(mesh, coeffsDict_)),
    rootNode_(0),
    terminalNodes_(),
    nodeLocations_(),
    terminalLocations_(),
    rootStartTime_(GREAT),
    rootDuration_(0.0),
    rootIntensity_(0.0),
    chi_(coeffsDict_.get<scalar>("chi")),
    Cm_(coeffsDict_.get<scalar>("cm")),
    Vm1D_(),
    Iion1D_(),
    activationTime_(),
    ionicModelPtr_(nullptr),
    verificationModelPtr_(nullptr),
    localStartNode_(-1),
    nLocalNodes_(0),
    terminalCurrent_(),
    terminalSource_(),
    outputPtr_(),
    exportVars_(),
    debugVars_(),
    reportSetup_(coeffsDict_.lookupOrDefault<Switch>("reportSetup", false)),
    pvdTimes_(),
    pvdFiles_()
{
    readGraphFile(coeffsDict_);
    readRootStimulus(coeffsDict_);
    initialiseState(initialDeltaT);
    initialiseOutputControls();

    if (coeffsDict_.found("verificationModel"))
    {
        verificationModelPtr_ = graphVerificationModel::New
        (
            coeffsDict_.subDict("verificationModel")
        );
    }
    openOutputFile();

    if (reportSetup_)
    {
        Info<< "conductionSystemDomain constructed as graph Purkinje model with "
            << graph_.nNodes << " nodes and " << graph_.nEdges << " edges."
            << nl << endl;
    }
}


void conductionSystemDomain::preProcess()
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->preProcess
        (
            time(),
            ionicModelPtr_.ptr(),
            Vm1D_,
            nodeLocations_
        );
    }
}


void conductionSystemDomain::advance(scalar t0, scalar dt)
{
    if (t0 == 0.0)
    {
        preProcess();
    }

    solverPtr_->advance(*this, t0, dt);
}


void conductionSystemDomain::assembleAppliedCurrent
(
    scalar t0,
    scalarField& appliedCurrent
) const
{
    appliedCurrent = 0.0;

    if (t0 >= rootStartTime_ && t0 <= (rootStartTime_ + rootDuration_))
    {
        appliedCurrent[rootNode_] += rootIntensity_;
    }

    forAll(terminalNodes_, i)
    {
        appliedCurrent[terminalNodes_[i]] -= terminalCurrent_[i];
    }
}


void conductionSystemDomain::reportAdvanceDiagnostics
(
    scalar t0,
    scalar dt
) const
{
    if (!debugVars_.size() || time().timeIndex() % 10 != 0)
    {
        return;
    }

    Info<< "[PurkinjeGraph] timeIndex=" << time().timeIndex()
        << " t=" << t0
        << " dt=" << dt;

    if (debugVars_.found("Vm"))
    {
        Info<< " Vm1D[min,max]=[" << gMin(Vm1D_) << ", " << gMax(Vm1D_) << "]";
    }

    if (debugVars_.found("Iion"))
    {
        Info<< " Iion1D[min,max]=[" << gMin(Iion1D_) << ", " << gMax(Iion1D_) << "]";
    }

    Info<< endl;
}


void conductionSystemDomain::terminalVm(scalarField& values) const
{
    values.setSize(terminalNodes_.size());
    values = 0.0;

    forAll(terminalNodes_, i)
    {
        values[i] = Vm1D_[terminalNodes_[i]];
    }
}


void conductionSystemDomain::terminalActivationTime(scalarField& values) const
{
    values.setSize(terminalNodes_.size());
    values = -1.0;

    forAll(terminalNodes_, i)
    {
        values[i] = activationTime_[terminalNodes_[i]];
    }
}


void conductionSystemDomain::setTerminalActivationTime(const scalarField& values)
{
    if (values.size() != terminalNodes_.size())
    {
        FatalErrorInFunction
            << "Expected " << terminalNodes_.size()
            << " terminal activation values but received "
            << values.size()
            << exit(FatalError);
    }

    forAll(terminalNodes_, i)
    {
        if (values[i] >= 0.0)
        {
            const label nodeI = terminalNodes_[i];
            if (activationTime_[nodeI] < 0.0 || values[i] < activationTime_[nodeI])
            {
                activationTime_[nodeI] = values[i];
            }
        }
    }
}


void conductionSystemDomain::setTerminalCoupling
(
    const scalarField& terminalCurrent,
    const scalarField& terminalSource
)
{
    if (terminalCurrent.size() != terminalNodes_.size())
    {
        FatalErrorInFunction
            << "Expected " << terminalNodes_.size()
            << " terminalCurrent values but received "
            << terminalCurrent.size()
            << exit(FatalError);
    }

    if (terminalSource.size() != terminalNodes_.size())
    {
        FatalErrorInFunction
            << "Expected " << terminalNodes_.size()
            << " terminalSource values but received "
            << terminalSource.size()
            << exit(FatalError);
    }

    terminalCurrent_ = terminalCurrent;
    terminalSource_ = terminalSource;
}


void conductionSystemDomain::end()
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->postProcess
        (
            time(),
            ionicModelPtr_.ptr(),
            Vm1D_,
            nodeLocations_
        );
    }
}


void conductionSystemDomain::write()
{
    if (!time().outputTime())
    {
        return;
    }

    DynamicList<scalar> values;
    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            forAll(Vm1D_, i)
            {
                values.append(Vm1D_[i]);
            }
        }
        else if (var == "Iion")
        {
            forAll(Iion1D_, i)
            {
                values.append(Iion1D_[i]);
            }
        }
        else if (var == "activationTime")
        {
            forAll(activationTime_, i)
            {
                values.append(activationTime_[i]);
            }
        }
        else if (var == "IcouplingSource")
        {
            forAll(terminalSource_, i)
            {
                values.append(terminalSource_[i]);
            }
        }
        else if (var == "IcouplingCurrent")
        {
            forAll(terminalCurrent_, i)
            {
                values.append(terminalCurrent_[i]);
            }
        }
    }

    if (Pstream::master() && outputPtr_.valid())
    {
        purkinjeModelIO::writeRow(outputPtr_.ref(), time().value(), values);
    }

    PtrList<scalarField> ionicFields(ionicExport_.size());

    if (ionicModelPtr_.valid() && ionicExport_.size() > 0)
    {
        const auto* statesPtr = ionicModelPtr_->ioStatesPtr();
        const auto* algebraicPtr = ionicModelPtr_->ioAlgebraicPtr();

        forAll(ionicExport_, i)
        {
            scalarField varField(graph_.nNodes, 0.0);
            label sIdx = ionicExportStateIndices_[i];
            label aIdx = ionicExportAlgebraicIndices_[i];

            if (sIdx >= 0 && statesPtr)
            {
                forAll(varField, nodeI)
                {
                    varField[nodeI] = (*statesPtr)[nodeI][sIdx];
                }
            }
            else if (aIdx >= 0 && algebraicPtr)
            {
                forAll(varField, nodeI)
                {
                    varField[nodeI] = (*algebraicPtr)[nodeI][aIdx];
                }
            }

            ionicFields.set(i, new scalarField(varField));

            if (Pstream::master() && ionicOutputPtrs_.set(i))
            {
                DynamicList<scalar> ionicValues(probeNodes_.size());
                forAll(probeNodes_, pI)
                {
                    ionicValues.append(varField[probeNodes_[pI]]);
                }
                purkinjeModelIO::writeRow(ionicOutputPtrs_[i], time().value(), ionicValues);
            }
        }
    }

    if (Pstream::master())
    {
        const fileName vtkDir
        (
            time().globalPath()/"postProcessing"/"purkinjeNetworkVTK"
        );

        char buf[32];
        snprintf(buf, sizeof(buf), "purkinjeNetwork_%06d.vtk", int(time().timeIndex()));
        const word vtkFilename(buf);

        wordList diagNames;
        PtrList<scalarField> diagFields;
        if (solverPtr_.valid())
        {
            solverPtr_->diagnosticFields(diagNames, diagFields);
        }

        const label nIonic = ionicExport_.size();
        const label nDiag = diagNames.size();

        wordList vtkNames(1 + nIonic + nDiag);
        PtrList<scalarField> vtkFields(1 + nIonic + nDiag);

        vtkNames[0] = "activationTime";
        vtkFields.set(0, new scalarField(activationTime_));

        forAll(ionicExport_, i)
        {
            vtkNames[1 + i] = ionicExport_[i];
            vtkFields.set(1 + i, new scalarField(ionicFields[i]));
        }

        forAll(diagNames, i)
        {
            vtkNames[1 + nIonic + i] = diagNames[i];
            vtkFields.set(1 + nIonic + i, new scalarField(diagFields[i]));
        }

        purkinjeModelIO::writeVTK
        (
            vtkDir,
            time().timeName(),
            time().timeIndex(),
            nodeLocations_,
            graph_.edgeNodeA,
            graph_.edgeNodeB,
            Vm1D_,
            Iion1D_,
            terminalNodes_,
            terminalSource_,
            vtkFields,
            vtkNames
        );

        pvdTimes_.append(time().value());
        pvdFiles_.append(vtkFilename);

        purkinjeModelIO::writeVTKSeries
        (
            vtkDir/"purkinjeNetwork.vtk.series",
            pvdTimes_,
            pvdFiles_
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
