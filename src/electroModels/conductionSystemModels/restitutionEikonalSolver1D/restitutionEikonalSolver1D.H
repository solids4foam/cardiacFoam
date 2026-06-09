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
    Foam::restitutionEikonalSolver1D

Description
    Dynamic, restitution-aware eikonal activation-time solver on the 1-D
    conduction graph. Unlike the static eikonalSolver1D, this solver keeps
    persistent recovery state and runs an event-driven loop each physical
    timestep: recover refractory nodes, schedule stimuli, propagate candidate
    arrival times from newly activated sources (refractory nodes block), then
    accept activations whose scheduled time falls in the step. APD and CV are
    drawn from a restitutionModel as functions of diastolic interval, so the
    solver reproduces rate-dependent conduction and conduction block, and
    reports where its assumptions break down (block / wavebreak / short DI).

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    restitutionEikonalSolver1D.C
\*---------------------------------------------------------------------------*/

#ifndef restitutionEikonalSolver1D_H
#define restitutionEikonalSolver1D_H

#include "conductionSystemSolver.H"
#include "restitutionModel.H"
#include "stimulusIO.H"
#include "fvMesh.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "labelList.H"
#include "autoPtr.H"
#include "Switch.H"

namespace Foam
{

class conductionSystemDomain;

/*---------------------------------------------------------------------------*\
                  Class restitutionEikonalSolver1D Declaration
\*---------------------------------------------------------------------------*/

class restitutionEikonalSolver1D
:
    public conductionSystemSolver
{
    enum nodeState { excitable = 0, refractory = 1 };

    autoPtr<restitutionModel> restitutionPtr_;

    // Persistent per-node recovery state.
    scalarField RT_;                 //- repolarization time
    scalarField DI_;                 //- diastolic interval at last activation
    scalarField APD_;                //- APD assigned at last activation
    scalarField nextTact_;           //- tentative next activation time
    labelList   state_;              //- excitable / refractory
    labelList   activatedFrom_;      //- neighbour that triggered activation (-1 = stimulus)
    labelList   nextTactSource_;     //- proposer of the current nextTact_ (-1 = stimulus)

    // Breakdown indicators (cumulative).
    labelList   blockCount_;
    labelList   wavebreakCount_;
    labelList   shortDICount_;
    scalarField minDI_;

    // Stimulus: shared S1-S2 protocol (timing) applied at node sites.
    labelList        stimSites_;
    StimulusProtocol stimProtocol_;

    scalar  initialAPD_;

    // Per-edge conduction heterogeneity: CV scaled by sqrt(g/gRef); an edge
    // with zero conductance is non-conducting (structural block, e.g. LBBB).
    Switch  useEdgeConductance_;
    scalar  refConductance_;

    Switch  reportSetup_;
    bool    initialised_;

    void initialiseState(const conductionSystemDomain& domain);

    restitutionEikonalSolver1D(const restitutionEikonalSolver1D&) = delete;
    void operator=(const restitutionEikonalSolver1D&) = delete;

public:

    OverrideTypeName("restitutionEikonalSolver1D");

    restitutionEikonalSolver1D
    (
        const fvMesh& mesh,
        const dictionary& solverCoeffs
    );

    virtual void advance
    (
        conductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) override;

    virtual void diagnosticFields
    (
        wordList& names,
        PtrList<scalarField>& fields
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
