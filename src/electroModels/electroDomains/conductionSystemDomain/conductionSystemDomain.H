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

Class
    Foam::ConductionSystemDomain

Description
    Concrete conduction-network domain (Purkinje network).

    Reads graph topology, PVJ mapping, root stimulus, and ionic model from
    the case dictionary.  Delegates time-stepping to a runtime-selectable
    conductionSystemSolver (monodomain1DSolver or eikonalSolver).

    Coupling exchange with the primary tissue domain is handled by a
    separate ElectroDomainCoupler (e.g. pvjResistanceCoupler).

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    conductionSystemDomain.C
\*---------------------------------------------------------------------------*/

#ifndef ConductionSystemDomain_H
#define ConductionSystemDomain_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "electroDomainInterface.H"
#include "electroDomainCouplingEndpoints.H"
#include "conductionGraph.H"
#include "conductionSystemSolver.H"
#include "ionicModel.H"
#include "purkinjeModelIO.H"

namespace Foam
{

class ConductionSystemDomain
:
    public electroDomainInterface,
    public networkCouplingEndpoint
{
    const fvMesh& supportMesh_;

    // ---- Graph topology (built from 'edges' list in dict) ----
    conductionGraph graph_;
    autoPtr<conductionSystemSolver> solverPtr_;

    // ---- PVJ mapping (parallel lists) ----
    labelList   pvjNodes_;        // Purkinje node indices that are PVJs
    pointField  pvjLocations_;    // 3D coordinates of each PVJ junction [m]
    //- Local cell ID of PVJ centroid (-1 if not on this proc)
    labelList   pvjCellIDs_;

    // ---- Root stimulus ----
    scalar      rootStartTime_;
    scalar      rootDuration_;
    scalar      rootIntensity_;   // [A/m^3] point-source applied to node 0

    // ---- Monodomain parameters for the Purkinje network ----
    scalar      chi_;             // surface-to-volume ratio [1/m]
    scalar      Cm_;              // membrane capacitance [F/m^2]

    // ---- Coupling ----
    scalar      pvjRadius_;       // sphere radius for current distribution [m]
    //- 3D cells within pvjRadius_ of each PVJ centroid
    List<labelList> pvjCellSets_;
    scalarField  pvjSphereVol_;   // total cell volume of each PVJ sphere [m^3]

    // ---- 1D state (size = graph_.nNodes) ----
    scalarField Vm1D_;            // transmembrane voltage [V]
    scalarField Iion1D_;          // ionic current [V/s]
    scalarField activationTime_;  // per-node activation time [s] (-1 = not yet activated)

    // ---- Ionic model (graph_.nNodes integration points) ----
    autoPtr<ionicModel> ionicModelPtr_;

    // Raw PVJ exchange term passed back from the coupling model.
    // Positive values remove charge from the network and inject it into tissue.
    // One entry per PVJ in the same order as pvjNodes_.
    scalarField terminalCurrent_;

    // Myocardium volumetric source equivalent associated with each PVJ.
    // Stored for diagnostics/output only, one entry per PVJ.
    scalarField terminalSource_;

    // ---- Output ----
    autoPtr<OFstream> outputPtr_;   // postProcessing/purkinjeNetwork.dat
    wordList          exportVars_;  // tokens written to purkinjeNetwork.dat
    wordList          debugVars_;   // tokens printed to terminal every 10 steps

    // ---- Private helpers ----
    void readPVJs(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);
    void openOutputFile();

    // No copy
    ConductionSystemDomain(const ConductionSystemDomain&) = delete;
    void operator=(const ConductionSystemDomain&) = delete;

public:

    TypeName("ConductionSystemDomain");

    ConductionSystemDomain
    (
        const fvMesh&     mesh,
        const dictionary& dict,
        scalar            initialDeltaT
    );

    virtual ~ConductionSystemDomain() = default;

    // ---- electroDomainInterface ----

    virtual const Time& time() const override
    {
        return supportMesh_.time();
    }

    virtual void advance(scalar t0, scalar dt) override;

    virtual void write() override;

    // ---- networkCouplingEndpoint ----

    virtual const labelList& terminalNodes() const override
    {
        return pvjNodes_;
    }

    virtual const pointField& terminalLocations() const override
    {
        return pvjLocations_;
    }

    virtual void terminalVm(scalarField& values) const override;

    virtual void setTerminalCoupling
    (
        const scalarField& terminalCurrent,
        const scalarField& terminalSource
    ) override;

    // ---- Graph accessors (used by conductionSystemSolver) ----

    scalarField& membranePotential() { return Vm1D_; }

    scalarField& ionicCurrent() { return Iion1D_; }

    ionicModel& ionicModelRef() { return ionicModelPtr_(); }

    const labelList& edgeStartNodes() const { return graph_.edgeNodeA; }

    const labelList& edgeEndNodes() const { return graph_.edgeNodeB; }

    const scalarField& edgeLengths() const { return graph_.edgeLengths; }

    const scalarField& edgeConductances() const { return graph_.edgeConductances; }

    scalar chi() const { return chi_; }

    scalar Cm() const { return Cm_; }

    const conductionGraph& graph() const { return graph_; }

    scalarField& activationTime() { return activationTime_; }

    void assembleAppliedCurrent
    (
        scalar t0,
        scalarField& appliedCurrent
    ) const;

    void reportAdvanceDiagnostics
    (
        scalar t0,
        scalar dt
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
