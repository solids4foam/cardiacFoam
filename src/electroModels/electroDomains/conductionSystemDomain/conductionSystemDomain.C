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

namespace Foam
{

defineTypeNameAndDebug(ConductionSystemDomain, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ConductionSystemDomain::readPVJs(const dictionary& dict)
{
    pvjNodes_     = labelList(dict.lookup("pvjNodes"));
    pvjLocations_ = pointField(dict.lookup("pvjLocations"));

    if (pvjNodes_.size() != pvjLocations_.size())
    {
        FatalErrorInFunction
            << "pvjNodes and pvjLocations must have the same size. "
            << "pvjNodes.size()=" << pvjNodes_.size()
            << " pvjLocations.size()=" << pvjLocations_.size()
            << exit(FatalError);
    }

    // Resolve each PVJ location to a local cell ID.
    // findCell returns -1 if the point is not in any local cell (parallel).
    pvjCellIDs_.setSize(pvjNodes_.size(), -1);
    forAll(pvjNodes_, k)
    {
        pvjCellIDs_[k] = supportMesh_.findCell(pvjLocations_[k]);
    }

    // Build pvjCellSets_: all LOCAL cells within pvjRadius_ of each PVJ
    // location.
    // pvjSphereVol_ accumulates per-processor; reduced globally below.
    const vectorField& CC = supportMesh_.C();
    pvjCellSets_.setSize(pvjNodes_.size());
    pvjSphereVol_.setSize(pvjNodes_.size(), 0.0);
    const scalarField& cellVols = supportMesh_.V();

    forAll(pvjNodes_, k)
    {
        const point& centre = pvjLocations_[k];
        DynamicList<label> cellsInRadius;
        forAll(CC, cellI)
        {
            if (mag(CC[cellI] - centre) <= pvjRadius_)
            {
                cellsInRadius.append(cellI);
                pvjSphereVol_[k] += cellVols[cellI];
            }
        }
        pvjCellSets_[k] = cellsInRadius;
    }

    // Sum sphere volumes across all processors to get global totals.
    forAll(pvjNodes_, k)
    {
        reduce(pvjSphereVol_[k], sumOp<scalar>());
    }

    Info<< "Purkinje PVJs: " << pvjNodes_.size() << " junctions, "
        << "pvjRadius=" << pvjRadius_*1e3 << " mm\n";
    forAll(pvjNodes_, k)
    {
        Info<< "  pvj" << k << " (node " << pvjNodes_[k]
            << ", location " << pvjLocations_[k]
            << ", localCentreCell=" << pvjCellIDs_[k]
            << "): " << pvjCellSets_[k].size()
            << " local cells, V_sphere_global="
            << pvjSphereVol_[k] << " m^3" << endl;
    }
}


void ConductionSystemDomain::openOutputFile()
{
    // Only master creates the file; workers leave outputPtr_ null.
    if (!Pstream::master()) return;

    // Build column headers from exportVars_.
    // Tokens: "Vm"        -> node0_Vm_mV  ... nodeN_Vm_mV
    //         "Iion"      -> node0_Iion   ... nodeN_Iion
    //         "Icoupling" / "IcouplingSource"
    //                      -> pvj0_IcouplingSource_Am3 ... pvjK_...
    //         "IcouplingCurrent"
    //                      -> pvj0_IcouplingCurrent ... pvjK_...
    DynamicList<word> colNames;
    colNames.append("time");
    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            for (label nI = 0; nI < graph_.nNodes; nI++)
                colNames.append("node" + Foam::name(nI) + "_Vm_mV");
        }
        else if (var == "Iion")
        {
            for (label nI = 0; nI < graph_.nNodes; nI++)
                colNames.append("node" + Foam::name(nI) + "_Iion");
        }
        else if (var == "Icoupling" || var == "IcouplingSource")
        {
            forAll(pvjNodes_, k)
            {
                colNames.append
                (
                    "pvj" + Foam::name(k) + "_IcouplingSource_Am3"
                );
            }
        }
        else if (var == "IcouplingCurrent")
        {
            forAll(pvjNodes_, k)
            {
                colNames.append
                (
                    "pvj" + Foam::name(k) + "_IcouplingCurrent"
                );
            }
        }
        else
        {
            WarningInFunction
                << "Unknown outputVariables token '" << var
                << "'. Valid tokens: Vm, Iion, Icoupling, "
                << "IcouplingSource, IcouplingCurrent." << endl;
        }
    }

    // Use the case root (one level up from processor*) so the file lands in
    // <case>/postProcessing/purkinjeNetwork.dat.
    const fileName outDir(time().path() / ".." / "postProcessing");
    outputPtr_ = purkinjeModelIO::openTimeSeries
    (
        outDir, "purkinjeNetwork.dat", colNames
    );

    Info<< "ConductionSystemDomain: writing to "
        << outDir / "purkinjeNetwork.dat" << endl;
}


void ConductionSystemDomain::readRootStimulus(const dictionary& dict)
{
    const dictionary& rsDict = dict.subDict("rootStimulus");
    rootStartTime_ = rsDict.get<scalar>("startTime");
    rootDuration_  = rsDict.get<scalar>("duration");
    rootIntensity_ = rsDict.get<scalar>("intensity");

    Info<< "Purkinje root stimulus: start=" << rootStartTime_
        << " duration=" << rootDuration_
        << " intensity=" << rootIntensity_ << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

ConductionSystemDomain::ConductionSystemDomain
(
    const fvMesh&     mesh,
    const dictionary& dict,
    const scalar      initialDeltaT
)
:
    supportMesh_(mesh),
    graph_(),
    solverPtr_(conductionSystemSolver::New(dict)),
    pvjNodes_(),
    pvjLocations_(),
    pvjCellIDs_(),
    rootStartTime_(0),
    rootDuration_(0),
    rootIntensity_(0),
    chi_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("chi")),
    Cm_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("cm")),
    pvjRadius_(dict.getOrDefault<scalar>("pvjRadius", 0.5e-3)),
    pvjCellSets_(),
    pvjSphereVol_(),
    Vm1D_(),
    Iion1D_(),
    activationTime_(),
    ionicModelPtr_(),
    terminalCurrent_(),
    terminalSource_(),
    exportVars_(),
    debugVars_()
{
    graph_.readFromDict(dict);
    readPVJs(dict);
    readRootStimulus(dict);

    // Initialise 1D membrane potential field.
    // Default -0.084 V matches the 3-D MyocardiumDomain.
    const scalar Vm1D_rest =
        dict.subDict("purkinjeNetworkModelCoeffs")
            .getOrDefault<scalar>("vm1DRest", -0.084);
    Vm1D_.setSize(graph_.nNodes, Vm1D_rest);
    Iion1D_.setSize(graph_.nNodes, 0.0);

    // Construct ionic model with graph_.nNodes integration points
    const dictionary& coeffsDict =
        dict.subDict("purkinjeNetworkModelCoeffs");

    ionicModelPtr_ = ionicModel::New
    (
        coeffsDict,
        graph_.nNodes,
        initialDeltaT,
        // solveVmWithinODESolver: always false for 1D PN
        // (Vm governed by diffusion PDE)
        false
    );

    Info<< "ConductionSystemDomain constructed with ionic model "
        << ionicModelPtr_->type() << endl;

    // Parse output variable selection (must precede openOutputFile)
    {
        const dictionary& ovDict =
            coeffsDict.subOrEmptyDict("outputVariables");
        exportVars_ = ovDict.getOrDefault<wordList>(
            "export", wordList{"Vm", "Icoupling"});
        debugVars_  = ovDict.getOrDefault<wordList>(
            "debug", wordList());
    }

    terminalCurrent_.setSize(pvjNodes_.size(), 0.0);
    terminalSource_.setSize(pvjNodes_.size(), 0.0);

    openOutputFile();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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
        appliedCurrent[0] += rootIntensity_;
    }

    forAll(pvjNodes_, k)
    {
        appliedCurrent[pvjNodes_[k]] -= terminalCurrent_[k];
    }
}


void ConductionSystemDomain::reportAdvanceDiagnostics
(
    scalar          t0,
    scalar          dt
) const
{
    // Diagnostic print — only if debugVars_ non-empty, every 10 timesteps.
    if (debugVars_.size() && time().timeIndex() % 10 == 0)
    {
        Info<< "[PK] timeIndex=" << time().timeIndex()
            << "  t=" << t0;
        if (debugVars_.found("Vm"))
        {
            Info<< "  Vm1D[mV]:";
            forAll(Vm1D_, nI)
                Info<< " n" << nI << "=" << Vm1D_[nI]*1e3;
        }
        if (debugVars_.found("Iion"))
        {
            Info<< "  Iion1D[0]=" << Iion1D_[0];
        }
        if
        (
            debugVars_.found("Icoupling")
         || debugVars_.found("IcouplingSource")
        )
        {
            Info<< "  IcouplingSource[Am3]:";
            forAll(pvjNodes_, k)
            {
                Info<< " pvj" << k << "=" << terminalSource_[k];
            }
        }
        if (debugVars_.found("IcouplingCurrent"))
        {
            Info<< "  IcouplingCurrent:";
            forAll(pvjNodes_, k)
            {
                Info<< " pvj" << k << "=" << terminalCurrent_[k];
            }
        }
        Info<< endl;
    }

    (void)dt;
}


void ConductionSystemDomain::terminalVm(scalarField& values) const
{
    values.setSize(pvjNodes_.size());
    values = 0.0;

    forAll(pvjNodes_, k)
    {
        values[k] = Vm1D_[pvjNodes_[k]];
    }
}


void ConductionSystemDomain::setTerminalCoupling
(
    const scalarField& terminalCurrent,
    const scalarField& terminalSource
)
{
    if (terminalCurrent.size() != pvjNodes_.size())
    {
        FatalErrorInFunction
            << "Expected " << pvjNodes_.size()
            << " terminalCurrent values but received "
            << terminalCurrent.size()
            << exit(FatalError);
    }

    if (terminalSource.size() != pvjNodes_.size())
    {
        FatalErrorInFunction
            << "Expected " << pvjNodes_.size()
            << " terminalSource values but received "
            << terminalSource.size()
            << exit(FatalError);
    }

    terminalCurrent_ = terminalCurrent;
    terminalSource_ = terminalSource;
}


void ConductionSystemDomain::write()
{
    if (!time().outputTime()) return;

    // Assemble selected columns — same token expansion as openOutputFile().
    DynamicList<scalar> vals;
    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            for (label nI = 0; nI < graph_.nNodes; nI++)
                vals.append(Vm1D_[nI] * 1e3);           // [mV]
        }
        else if (var == "Iion")
        {
            for (label nI = 0; nI < graph_.nNodes; nI++)
                vals.append(Iion1D_[nI]);
        }
        else if (var == "Icoupling" || var == "IcouplingSource")
        {
            forAll(pvjNodes_, k)
                vals.append(terminalSource_[k]);
        }
        else if (var == "IcouplingCurrent")
        {
            forAll(pvjNodes_, k)
                vals.append(terminalCurrent_[k]);
        }
    }

    if (Pstream::master())
    {
        purkinjeModelIO::writeRow(outputPtr_.ref(), time().value(), vals);
    }
}


} // End namespace Foam

// ************************************************************************* //
