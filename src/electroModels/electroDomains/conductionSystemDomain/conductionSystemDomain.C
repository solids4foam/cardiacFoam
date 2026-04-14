/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "conductionSystemDomain.H"

#include "DynamicList.H"
#include "IOdictionary.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(ConductionSystemDomain, 0);

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

    return solverType != "eikonalSolver1D" && solverType != "eikonalSolver";
}

} // End anonymous namespace


autoPtr<ConductionSystemDomain> ConductionSystemDomain::New
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

    if
    (
        domainType != "purkinjeGraphModel"
     && domainType != "conductionSystemDomain"
    )
    {
        FatalErrorInFunction
            << "Unsupported conductionSystemDomain type '" << domainType
            << "'. Supported graph domain types are 'purkinjeGraphModel' and "
            << "'conductionSystemDomain'."
            << exit(FatalError);
    }

    return autoPtr<ConductionSystemDomain>
    (
        new ConductionSystemDomain(mesh, dict, initialDeltaT)
    );
}


void ConductionSystemDomain::readGraphFile(const dictionary& dict)
{
    // graphFile is optional; if not present, skip graph reading (no Purkinje network)
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

    graph_.readFromDict(graphDict);

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


void ConductionSystemDomain::readRootStimulus(const dictionary& dict)
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


void ConductionSystemDomain::initialiseState(const scalar initialDeltaT)
{
    if (selectedSolverRequiresIonicModel(coeffsDict_))
    {
        ionicModelPtr_ = ionicModel::New
        (
            coeffsDict_,
            graph_.nNodes,
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
        Info<< "ConductionSystemDomain ionic model: "
            << ionicModelPtr_->type() << nl << endl;
    }
}


void ConductionSystemDomain::initialiseOutputControls()
{
    const dictionary& ovDict = coeffsDict_.subOrEmptyDict("outputVariables");

    exportVars_ = ovDict.getOrDefault<wordList>
    (
        "export",
        wordList{"Vm", "Icoupling"}
    );

    debugVars_ = ovDict.getOrDefault<wordList>
    (
        "debug",
        wordList()
    );
}


void ConductionSystemDomain::openOutputFile()
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
        else if (var == "Icoupling" || var == "IcouplingSource")
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

    Info<< "ConductionSystemDomain: writing to "
        << outDir/"purkinjeNetwork.dat" << nl << endl;
}


ConductionSystemDomain::ConductionSystemDomain
(
    const fvMesh& mesh,
    const dictionary& dict,
    const scalar initialDeltaT
)
:
    supportMesh_(mesh),
    coeffsDict_
    (
        dict.found("purkinjeGraphModelCoeffs")
      ? dict.subDict("purkinjeGraphModelCoeffs")
      : dict.found("purkinjeNetworkModelCoeffs")
      ? dict.subDict("purkinjeNetworkModelCoeffs")
      : dict
    ),
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
    ionicModelPtr_(),
    terminalCurrent_(),
    terminalSource_(),
    outputPtr_(),
    exportVars_(),
    debugVars_(),
    pvdTimes_(),
    pvdFiles_()
{
    readGraphFile(coeffsDict_);
    readRootStimulus(coeffsDict_);
    initialiseState(initialDeltaT);
    initialiseOutputControls();
    openOutputFile();

    Info<< "ConductionSystemDomain constructed as graph Purkinje model with "
        << graph_.nNodes << " nodes and " << graph_.nEdges << " edges."
        << nl << endl;
}


void ConductionSystemDomain::advance
(
    scalar t0,
    scalar dt
)
{
    solverPtr_->advance(*this, t0, dt);
}


void ConductionSystemDomain::assembleAppliedCurrent
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


void ConductionSystemDomain::reportAdvanceDiagnostics
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


void ConductionSystemDomain::terminalVm(scalarField& values) const
{
    values.setSize(terminalNodes_.size());
    values = 0.0;

    forAll(terminalNodes_, i)
    {
        values[i] = Vm1D_[terminalNodes_[i]];
    }
}


void ConductionSystemDomain::terminalActivationTime(scalarField& values) const
{
    values.setSize(terminalNodes_.size());
    values = -1.0;

    forAll(terminalNodes_, i)
    {
        values[i] = activationTime_[terminalNodes_[i]];
    }
}


void ConductionSystemDomain::setTerminalCoupling
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


void ConductionSystemDomain::write()
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
        else if (var == "Icoupling" || var == "IcouplingSource")
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

    if (Pstream::master())
    {
        const fileName vtkDir
        (
            time().globalPath()/"postProcessing"/"purkinjeNetworkVTK"
        );

        char buf[32];
        snprintf(buf, sizeof(buf), "purkinjeNetwork_%06d.vtk", int(time().timeIndex()));
        const word vtkFilename(buf);

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
            terminalSource_
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
