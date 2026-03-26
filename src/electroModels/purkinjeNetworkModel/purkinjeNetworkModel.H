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
    labelList   pvjCellIDs_;      // matching 3D mesh cell indices

    // ---- Root stimulus ----
    scalar      rootStartTime_;
    scalar      rootDuration_;
    scalar      rootIntensity_;   // [A/m^3] point-source applied to node 0

    // ---- Monodomain parameters for the Purkinje network ----
    scalar      chi_;             // surface-to-volume ratio [1/m]
    scalar      Cm_;              // membrane capacitance [F/m^2]

    // ---- Coupling ----
    scalar      R_pvj_;           // global PVJ resistance [Ohm.m^2]

    // ---- 1D state (size = nNodes_) ----
    scalarField Vm1D_;            // transmembrane voltage [V]
    scalarField Iion1D_;          // ionic current [V/s]

    // ---- Ionic model (nNodes_ integration points) ----
    autoPtr<ionicModel> ionicModelPtr_;

    // ---- Private helpers ----
    void readGraph(const dictionary& dict);
    void readPVJs(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
