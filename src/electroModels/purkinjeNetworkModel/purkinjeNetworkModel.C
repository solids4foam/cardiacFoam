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
    pvjNodes_   = labelList(dict.lookup("pvjNodes"));
    pvjCellIDs_ = labelList(dict.lookup("pvjCellIDs"));

    if (pvjNodes_.size() != pvjCellIDs_.size())
    {
        FatalErrorInFunction
            << "pvjNodes and pvjCellIDs must have the same size. "
            << "pvjNodes.size()=" << pvjNodes_.size()
            << " pvjCellIDs.size()=" << pvjCellIDs_.size()
            << exit(FatalError);
    }

    const label nCells = mesh_.nCells();
    forAll(pvjCellIDs_, k)
    {
        if (pvjCellIDs_[k] < 0 || pvjCellIDs_[k] >= nCells)
        {
            FatalErrorInFunction
                << "pvjCellIDs[" << k << "] = " << pvjCellIDs_[k]
                << " is out of range [0, " << nCells - 1 << "]."
                << exit(FatalError);
        }
    }

    Info<< "Purkinje PVJs: " << pvjNodes_.size() << " junctions." << endl;
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
    pvjCellIDs_(),
    rootStartTime_(0),
    rootDuration_(0),
    rootIntensity_(0),
    chi_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("chi")),
    Cm_(dict.subDict("purkinjeNetworkModelCoeffs").get<scalar>("Cm")),
    R_pvj_(dict.get<scalar>("R_pvj")),
    Vm1D_(),
    Iion1D_(),
    ionicModelPtr_()
{
    readGraph(dict);
    readPVJs(dict);
    readRootStimulus(dict);

    // Initialise 1D state
    Vm1D_.setSize(nNodes_, -0.084);   // Stewart resting potential [V]
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
        << ionicModelPtr_->type() << nl << endl;
}


// * * * * * * * * * * * * evolve (stub for Task 3) * * * * * * * * * * * * //

void purkinjeNetworkModel::evolve
(
    scalar          /*t0*/,
    scalar          /*dt*/,
    volScalarField& /*externalStimulusCurrent*/
)
{
    // Implemented in Task 3
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
