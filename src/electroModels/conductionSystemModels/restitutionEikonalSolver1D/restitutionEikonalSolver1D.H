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
    conduction graph. Keeps persistent recovery state and runs an event-driven
    loop each physical timestep: recover refractory nodes, schedule stimuli,
    propagate candidate arrival times from newly activated sources (refractory
    nodes block), then accept activations whose scheduled time falls in the
    step. APD and CV are drawn from a restitutionModel as functions of
    diastolic interval. Reports breakdown indicators (block / wavebreak /
    short DI). Records multi-beat history via activationHistory.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    restitutionEikonalSolver1D.C
\*---------------------------------------------------------------------------*/

#ifndef restitutionEikonalSolver1D_H
#define restitutionEikonalSolver1D_H

#include "conductionSystemSolver.H"
#include "restitutionModel.H"
#include "activationHistory.H"
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

    scalarField RT_;
    scalarField DI_;
    scalarField APD_;
    scalarField nextTact_;
    labelList   state_;
    labelList   activatedFrom_;
    labelList   nextTactSource_;

    labelList   blockCount_;
    labelList   wavebreakCount_;
    labelList   shortDICount_;
    scalarField minDI_;

    activationHistory history_;

    labelList        stimSites_;
    StimulusProtocol stimProtocol_;

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
