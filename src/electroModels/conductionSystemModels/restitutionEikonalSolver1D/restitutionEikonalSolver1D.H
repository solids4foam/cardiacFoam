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
    Restitution-aware eikonal activation-time solver for the 1D conduction
    graph. State is tracked via lastActTime per node; all regimes are expressed
    as beat-to-beat intervals (beatInterval = te - lastActTime):
      beatInterval < minBeatInterval   — refractory / relative refractory (block)
      beatInterval >= minBeatInterval  — excitable, activated by incoming wave
      beatInterval >= escapeInterval   — spontaneous firing (funny current)
    minBeatInterval = apdNominal + minDI. CV is looked up as a function of
    DI = beatInterval - apdNominal. Reports breakdown indicators and records
    multi-beat history via activationHistory.

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
    autoPtr<restitutionModel> restitutionPtr_;

    scalarField lastActTime_;   //- time of last activation [s]; -GREAT if never fired
    scalarField DI_;
    scalarField nextTact_;
    labelList   activatedFrom_;
    labelList   nextTactSource_;

    labelList   blockCount_;
    labelList   wavebreakCount_;
    scalarField minDI_;

    activationHistory history_;

    labelList        stimSites_;
    StimulusProtocol stimProtocol_;

    scalar  apdNominal_;        //- constant APD [s] (dict: apdNominal)
    scalar  minBeatInterval_;   //- apdNominal + minDI: minimum beat-to-beat interval [s]
    scalar  escapeInterval_;    //- spontaneous firing beat-to-beat interval [s] (dict: escapeInterval)
    scalar  tStart_;            //- simulation start time [s]

    Switch  useEdgeConductance_;
    scalar  refConductance_;

    Switch  reportSetup_;
    bool    initialised_;

    void initialiseState(const conductionSystemDomain& domain, scalar t0);

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
