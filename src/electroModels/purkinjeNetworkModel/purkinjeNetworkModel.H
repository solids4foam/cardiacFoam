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

#include "purkinjeModel.H"
#include "ionicModel.H"
#include "ecgModelIO.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class purkinjeNetworkModel Declaration
\*---------------------------------------------------------------------------*/

class purkinjeNetworkModel
:
    public purkinjeModel
{
    // ---- Graph topology (built from 'edges' list in dict) ----
    label       nNodes_;          // total node count in the Purkinje tree
    label       nEdges_;
    labelList   edgeNodeA_;       // first endpoint of each edge
    labelList   edgeNodeB_;       // second endpoint of each edge
    scalarField edgeLength_;      // L  [m]
    scalarField edgeConductance_; // sigma [S/m]

    // ---- PVJ mapping (parallel lists) ----
    labelList   pvjNodes_;        // Purkinje node indices that are PVJs
    pointField  pvjLocations_;    // 3D coordinates of each PVJ junction [m]
    labelList   pvjCellIDs_;      // local cell ID of PVJ centroid (-1 if not on this proc)

    // ---- Root stimulus ----
    scalar      rootStartTime_;
    scalar      rootDuration_;
    scalar      rootIntensity_;   // [A/m^3] point-source applied to node 0

    // ---- Monodomain parameters for the Purkinje network ----
    scalar      chi_;             // surface-to-volume ratio [1/m]
    scalar      Cm_;              // membrane capacitance [F/m^2]

    // ---- Coupling ----
    scalar      R_pvj_;           // global PVJ resistance [Ohm.m^2]
    scalar          pvjRadius_;       // sphere radius for current distribution [m]
    List<labelList> pvjCellSets_;    // 3D cells within pvjRadius_ of each PVJ centroid
    scalarField     pvjSphereVol_;   // total cell volume of each PVJ sphere [m^3]

    // ---- 1D state (size = nNodes_) ----
    scalarField Vm1D_;            // transmembrane voltage [V]
    // Ionic source term [V/s]: stores Iion/(chi*Cm) — same scaling convention
    // as Iion_ in monoDomainElectro (ionicModel::solveODE output is already scaled).
    scalarField Iion1D_;

    // ---- Ionic model (nNodes_ integration points) ----
    autoPtr<ionicModel> ionicModelPtr_;

    // ---- Output (mirrors ecgModelIO pattern) ----
    autoPtr<OFstream> outputPtr_;   // postProcessing/purkinjeNetwork.dat
    wordList          exportVars_;  // tokens written to purkinjeNetwork.dat
    wordList          debugVars_;   // tokens printed to terminal every 10 steps

    // ---- Private helpers ----
    void readGraph(const dictionary& dict);
    void readPVJs(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);
    void openOutputFile();

    // No copy
    purkinjeNetworkModel(const purkinjeNetworkModel&) = delete;
    void operator=(const purkinjeNetworkModel&) = delete;

public:

    TypeName("purkinjeNetworkModel");

    // Constructors

        purkinjeNetworkModel
        (
            const volScalarField& Vm,
            const dictionary&     dict,
            const scalar          initialDeltaT
        );

    //- Destructor
    virtual ~purkinjeNetworkModel() = default;

    // Member Functions

        virtual void evolve
        (
            scalar          t0,
            scalar          dt,
            volScalarField& externalStimulusCurrent
        ) override;

        virtual void write() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
