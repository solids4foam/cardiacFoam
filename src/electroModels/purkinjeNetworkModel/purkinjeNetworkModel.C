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

#include "purkinjeNetworkModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(purkinjeNetworkModel, 0);
addToRunTimeSelectionTable(purkinjeModel, purkinjeNetworkModel, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void purkinjeNetworkModel::readGraph(const dictionary& dict)
{
    // Each entry in 'edges' is a list of 4 scalars:
    // ( nodeA  nodeB  length  conductance )
    const List<scalarList> edgeEntries(dict.get<List<scalarList>>("edges"));
    nEdges_ = edgeEntries.size();

    edgeNodeA_.setSize(nEdges_);
    edgeNodeB_.setSize(nEdges_);
    edgeLength_.setSize(nEdges_);
    edgeConductance_.setSize(nEdges_);

    label maxNode = 0;
    forAll(edgeEntries, eI)
    {
        const scalarList& e = edgeEntries[eI];
        if (e.size() != 4)
        {
            FatalErrorInFunction
                << "Each entry in 'edges' must have 4 values: "
                << "(nodeA nodeB length conductance). Entry " << eI
                << " has " << e.size() << " values."
                << exit(FatalError);
        }
        edgeNodeA_[eI]       = static_cast<label>(e[0]);
        edgeNodeB_[eI]       = static_cast<label>(e[1]);
        edgeLength_[eI]      = e[2];
        edgeConductance_[eI] = e[3];
        maxNode = max(maxNode, max(edgeNodeA_[eI], edgeNodeB_[eI]));
    }

    nNodes_ = maxNode + 1;

    Info<< "Purkinje network: " << nNodes_ << " nodes, "
        << nEdges_ << " edges." << endl;
}

void purkinjeNetworkModel::readPVJs(const dictionary& dict)
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
        pvjCellIDs_[k] = mesh_.findCell(pvjLocations_[k]);
    }

    // Build pvjCellSets_: all LOCAL cells within pvjRadius_ of each PVJ location.
    // pvjSphereVol_ accumulates per-processor; reduced globally below.
    const vectorField& CC = mesh_.C();
    pvjCellSets_.setSize(pvjNodes_.size());
    pvjSphereVol_.setSize(pvjNodes_.size(), 0.0);
    const scalarField& cellVols = mesh_.V();

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
            << " local cells, V_sphere_global=" << pvjSphereVol_[k] << " m^3" << endl;
    }
}

void purkinjeNetworkModel::openOutputFile()
{
    // Only master creates the file; workers leave outputPtr_ null.
    if (!Pstream::master()) return;

    // Build column headers from exportVars_.
    // Tokens: "Vm"        -> node0_Vm_mV  ... nodeN_Vm_mV
    //         "Iion"      -> node0_Iion   ... nodeN_Iion
    //         "Icoupling" -> pvj0_Icoupling_Am3 ... pvjK_Icoupling_Am3
    DynamicList<word> colNames;
    colNames.append("time");
    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            for (label nI = 0; nI < nNodes_; nI++)
                colNames.append("node" + Foam::name(nI) + "_Vm_mV");
        }
        else if (var == "Iion")
        {
            for (label nI = 0; nI < nNodes_; nI++)
                colNames.append("node" + Foam::name(nI) + "_Iion");
        }
        else if (var == "Icoupling")
        {
            forAll(pvjNodes_, k)
                colNames.append("pvj" + Foam::name(k) + "_Icoupling_Am3");
        }
        else
        {
            WarningInFunction
                << "Unknown outputVariables token '" << var
                << "'. Valid tokens: Vm, Iion, Icoupling." << endl;
        }
    }

    // Use the case root (one level up from processor*) so the file lands in
    // <case>/postProcessing/purkinjeNetwork.dat.
    const fileName outDir(mesh_.time().path() / ".." / "postProcessing");
    outputPtr_ = ecgModelIO::openTimeSeries(outDir, "purkinjeNetwork.dat", colNames);

    Info<< "purkinjeNetworkModel: writing to "
        << outDir / "purkinjeNetwork.dat" << endl;
}

void purkinjeNetworkModel::readRootStimulus(const dictionary& dict)
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

purkinjeNetworkModel::purkinjeNetworkModel
(
    const volScalarField& Vm,
    const dictionary&     dict,
    const scalar          initialDeltaT
)
:
    purkinjeModel(Vm, dict),
    nNodes_(0),
    nEdges_(0),
    edgeNodeA_(),
    edgeNodeB_(),
    edgeLength_(),
    edgeConductance_(),
    pvjNodes_(),
    pvjLocations_(),
    pvjCellIDs_(),
    rootStartTime_(0),
    rootDuration_(0),
    rootIntensity_(0),
    chi_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("chi")),
    Cm_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("Cm")),
    R_pvj_(dict.get<scalar>("R_pvj")),
    pvjRadius_(dict.getOrDefault<scalar>("pvjRadius", 0.5e-3)),
    pvjCellSets_(),
    pvjSphereVol_(),
    Vm1D_(),
    Iion1D_(),
    ionicModelPtr_(),
    exportVars_(),
    debugVars_()
{
    readGraph(dict);
    readPVJs(dict);
    readRootStimulus(dict);

    // Initialise 1D state.
    // Resting value depends on the ionic model:
    //   Stewart/TNNP/etc: physical voltage ~-0.084 V
    //   BuenoOrovio:      normalised u = 0.0 (dimensionless, rest = 0)
    // Read from dict with a default that works for dimensionless models.
    const scalar Vm1D_rest =
        dict.subDict("purkinjeNetworkModelCoeffs")
            .getOrDefault<scalar>("Vm1D_rest", 0.0);
    Vm1D_.setSize(nNodes_, Vm1D_rest);
    Iion1D_.setSize(nNodes_, 0.0);

    // Construct ionic model with nNodes_ integration points
    const dictionary& coeffsDict =
        dict.subDict("purkinjeNetworkModelCoeffs");

    ionicModelPtr_ = ionicModel::New
    (
        coeffsDict,
        nNodes_,
        initialDeltaT,
        false   // solveVmWithinODESolver: always false for 1D PN (Vm governed by diffusion PDE)
    );

    Info<< "purkinjeNetworkModel constructed with ionic model "
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

    openOutputFile();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void purkinjeNetworkModel::evolve
(
    scalar          t0,
    scalar          dt,
    volScalarField& externalStimulusCurrent
)
{
    // --- Step 1: Apply root stimulus at node 0 ---
    scalarField rootCurrent(nNodes_, 0.0);
    if (t0 >= rootStartTime_ && t0 <= (rootStartTime_ + rootDuration_))
    {
        rootCurrent[0] = rootIntensity_;
    }

    // --- Step 2: Advance ionic model ---
    ionicModelPtr_->solveODE(t0, dt, Vm1D_, Iion1D_);

    // --- Step 3: Graph Laplacian + explicit Euler ---
    scalarField laplacian(nNodes_, 0.0);
    for (label eI = 0; eI < nEdges_; eI++)
    {
        const label  A     = edgeNodeA_[eI];
        const label  B     = edgeNodeB_[eI];
        const scalar L     = edgeLength_[eI];
        const scalar sigma = edgeConductance_[eI];
        const scalar flux  = sigma * (Vm1D_[B] - Vm1D_[A]) / (L * L);
        laplacian[A] += flux;
        laplacian[B] -= flux;
    }

    const scalar dtOverChiCm = dt / (chi_ * Cm_);
    for (label nI = 0; nI < nNodes_; nI++)
    {
        Vm1D_[nI] += dtOverChiCm * (laplacian[nI] + rootCurrent[nI])
                   - dt * Iion1D_[nI];
    }

    // --- Step 4: Scatter coupling current to 3D externalStimulusCurrent ---
    // Vergara-Quarteroni bidirectional coupling:
    //   I_total [A/m^2] = (Vm1D_[pvj] - Vm3D_centre) / R_pvj_
    // To inject as volumetric current density [A/m^3] uniformly over the sphere:
    //   I_density [A/m^3] = I_total / V_sphere
    // This preserves total current and produces the correct magnitude to
    // activate the 3D domain (equivalent to monodomainStimulus intensity).
    scalarField& extI = externalStimulusCurrent.primitiveFieldRef();

    forAll(pvjNodes_, k)
    {
        // Get Vm3D at the PVJ centre. Only the processor that owns the cell
        // has a valid local ID; others contribute 0. maxOp gives the real value.
        scalar Vm3D_pvj = 0.0;
        if (pvjCellIDs_[k] >= 0)
        {
            Vm3D_pvj = Vm_[pvjCellIDs_[k]];
        }
        reduce(Vm3D_pvj, maxOp<scalar>());

        const scalar Icoupling =
            (Vm1D_[pvjNodes_[k]] - Vm3D_pvj)
            / (R_pvj_ * pvjSphereVol_[k]);

        forAll(pvjCellSets_[k], cI)
        {
            extI[pvjCellSets_[k][cI]] += Icoupling;
        }
    }

    externalStimulusCurrent.correctBoundaryConditions();

    // Diagnostic print — only if debugVars_ non-empty, every 10 timesteps.
    // Tokens: "Vm" -> per-node voltage, "Iion" -> Iion1D[0], "Icoupling" -> rootActive flag.
    if (debugVars_.size() && mesh_.time().timeIndex() % 10 == 0)
    {
        Info<< "[PK] timeIndex=" << mesh_.time().timeIndex()
            << "  t=" << t0;
        if (findIndex(debugVars_, word("Vm")) >= 0)
        {
            Info<< "  Vm1D[mV]:";
            forAll(Vm1D_, nI)
                Info<< " n" << nI << "=" << Vm1D_[nI]*1e3;
        }
        if (findIndex(debugVars_, word("Iion")) >= 0)
        {
            Info<< "  Iion1D[0]=" << Iion1D_[0];
        }
        if (findIndex(debugVars_, word("Icoupling")) >= 0)
        {
            Info<< "  rootActive="
                << label(t0 >= rootStartTime_ && t0 <= rootStartTime_ + rootDuration_);
        }
        Info<< endl;
    }
}


void purkinjeNetworkModel::write()
{
    if (!mesh_.time().outputTime()) return;

    // Gather Vm3D at all PVJ centres once (collective — all procs must participate).
    // maxOp: the owning processor provides the real value; others contribute 0.
    scalarField Vm3D_pvjs(pvjNodes_.size(), 0.0);
    forAll(pvjNodes_, k)
    {
        if (pvjCellIDs_[k] >= 0)
            Vm3D_pvjs[k] = Vm_[pvjCellIDs_[k]];
        reduce(Vm3D_pvjs[k], maxOp<scalar>());
    }

    // Assemble selected columns — same token expansion as openOutputFile().
    DynamicList<scalar> vals;
    for (const word& var : exportVars_)
    {
        if (var == "Vm")
        {
            for (label nI = 0; nI < nNodes_; nI++)
                vals.append(Vm1D_[nI] * 1e3);           // [mV]
        }
        else if (var == "Iion")
        {
            for (label nI = 0; nI < nNodes_; nI++)
                vals.append(Iion1D_[nI]);
        }
        else if (var == "Icoupling")
        {
            forAll(pvjNodes_, k)
                vals.append(
                    (Vm1D_[pvjNodes_[k]] - Vm3D_pvjs[k])
                    / (R_pvj_ * pvjSphereVol_[k])        // [A/m^3]
                );
        }
    }

    if (Pstream::master())
    {
        ecgModelIO::writeRow(outputPtr_.ref(), mesh_.time().value(), vals);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
