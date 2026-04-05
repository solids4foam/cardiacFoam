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
    Foam::purkinjeNetworkModel

Description
    Concrete Purkinje network model: reads graph topology, PVJ mapping,
    root stimulus and ionic model from the dictionary.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    purkinjeNetworkModel.C
\*---------------------------------------------------------------------------*/

#ifndef purkinjeNetworkModel_H
#define purkinjeNetworkModel_H

#include "graphConductionSystemDomain.H"
#include "conductionGraph.H"
#include "ionicModel.H"
#include "ecgModelIO.H"
#include "OFstream.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class purkinjeNetworkModel Declaration
\*---------------------------------------------------------------------------*/

class purkinjeNetworkModel
:
    public GraphConductionSystemDomain
{
    const fvMesh& supportMesh_;

    // ---- Graph topology (built from 'edges' list in dict) ----
    conductionGraph graph_;

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
    //- Effective PVJ resistance coefficient used in V/R coupling
    scalar      R_pvj_;
    scalar      pvjRadius_;       // sphere radius for current distribution [m]
    //- 3D cells within pvjRadius_ of each PVJ centroid
    List<labelList> pvjCellSets_;
    scalarField  pvjSphereVol_;   // total cell volume of each PVJ sphere [m^3]

    // ---- 1D state (size = graph_.nNodes) ----
    scalarField Vm1D_;            // transmembrane voltage [V]
    scalarField Iion1D_;           // ionic current [V/s]
    scalarField activationTime_;   // per-node activation time [s] (-1 = not yet activated)

    // ---- Ionic model (graph_.nNodes integration points) ----
    autoPtr<ionicModel> ionicModelPtr_;
    // Raw PVJ exchange term passed back from the coupling model.
    // Positive values remove charge from the network and inject it into tissue.
    // One entry per PVJ in the same order as pvjNodes_.
    scalarField terminalCurrent_;

    // Myocardium volumetric source equivalent associated with each PVJ.
    // Stored for diagnostics/output only, one entry per PVJ.
    scalarField terminalSource_;

    // ---- Output (mirrors ecgModelIO pattern) ----
    autoPtr<OFstream> outputPtr_;   // postProcessing/purkinjeNetwork.dat
    wordList          exportVars_;  // tokens written to purkinjeNetwork.dat
    wordList          debugVars_;   // tokens printed to terminal every 10 steps

    // ---- Private helpers ----
    void readPVJs(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);
    void openOutputFile();

    // No copy
    purkinjeNetworkModel(const purkinjeNetworkModel&) = delete;
    void operator=(const purkinjeNetworkModel&) = delete;

public:

    OverrideTypeName("purkinjeNetworkModel");

    // Constructors

        purkinjeNetworkModel
        (
            const ConductionSystemDomainContext& context,
            const dictionary&                    dict,
            const scalar                         initialDeltaT
        );

    //- Destructor
    virtual ~purkinjeNetworkModel() = default;

    // Member Functions

        virtual scalarField& membranePotential() override
        {
            return Vm1D_;
        }

        virtual scalarField& ionicCurrent() override
        {
            return Iion1D_;
        }

        virtual ionicModel& ionicModelRef() override
        {
            return ionicModelPtr_();
        }

        virtual const labelList& edgeStartNodes() const override
        {
            return graph_.edgeNodeA;
        }

        virtual const labelList& edgeEndNodes() const override
        {
            return graph_.edgeNodeB;
        }

        virtual const scalarField& edgeLengths() const override
        {
            return graph_.edgeLengths;
        }

        virtual const scalarField& edgeConductances() const override
        {
            return graph_.edgeConductances;
        }

        virtual scalar chi() const override
        {
            return chi_;
        }

        virtual scalar Cm() const override
        {
            return Cm_;
        }

        virtual const conductionGraph& graph() const override
        {
            return graph_;
        }

        virtual scalarField& activationTime() override
        {
            return activationTime_;
        }

        virtual void assembleAppliedCurrent
        (
            scalar t0,
            scalarField& appliedCurrent
        ) const override;

        virtual void reportAdvanceDiagnostics
        (
            scalar t0,
            scalar dt
        ) const override;

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

        virtual void write() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
