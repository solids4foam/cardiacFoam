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
    conductionSystemDomain

Description
    Graph-backed Purkinje/conduction-system domain.

SourceFiles
    conductionSystemDomain.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef conductionSystemDomain_H
#define conductionSystemDomain_H

#include "autoPtr.H"
#include "conductionGraph.H"
#include "conductionSystemSolver.H"
#include "electroDomainCouplingEndpoints.H"
#include "electroDomainInterface.H"
#include "fvMesh.H"
#include "ionicModel.H"
#include "labelList.H"
#include "OFstream.H"
#include "pointField.H"
#include "purkinjeModelIO.H"
#include "scalarField.H"
#include "Switch.H"

namespace Foam
{

class dictionary;

class conductionSystemDomain
:
    public electroDomainInterface,
    public networkCouplingEndpoint
{
    const fvMesh& supportMesh_;
    const dictionary& coeffsDict_;

    conductionGraph graph_;
    autoPtr<conductionSystemSolver> solverPtr_;

    label rootNode_;

    labelList terminalNodes_;
    pointField nodeLocations_;
    pointField terminalLocations_;

    scalar rootStartTime_;
    scalar rootDuration_;
    scalar rootIntensity_;

    scalar chi_;
    scalar Cm_;

    scalarField Vm1D_;
    scalarField Iion1D_;
    scalarField activationTime_;

    autoPtr<ionicModel> ionicModelPtr_;

    scalarField terminalCurrent_;
    scalarField terminalSource_;

    autoPtr<OFstream> outputPtr_;
    wordList exportVars_;
    wordList debugVars_;
    Switch reportSetup_;

    // PVD time-series index: parallel lists of physical time and VTK filename
    scalarList pvdTimes_;
    wordList   pvdFiles_;

    void readGraphFile(const dictionary& dict);
    void readRootStimulus(const dictionary& dict);
    void initialiseState(scalar initialDeltaT);
    void initialiseOutputControls();
    void openOutputFile();

    conductionSystemDomain(const conductionSystemDomain&) = delete;
    void operator=(const conductionSystemDomain&) = delete;

public:

    TypeName("conductionSystemDomain");

    static autoPtr<conductionSystemDomain> New
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar initialDeltaT
    );

    conductionSystemDomain
    (
        const fvMesh& mesh,
        const dictionary& dict,
        scalar initialDeltaT
    );

    virtual ~conductionSystemDomain() = default;

    virtual const Time& time() const override
    {
        return supportMesh_.time();
    }

    virtual void advance(scalar t0, scalar dt) override;

    virtual void write() override;

    virtual const labelList& terminalNodes() const override
    {
        return terminalNodes_;
    }

    virtual const pointField& terminalLocations() const override
    {
        return terminalLocations_;
    }

    virtual void terminalVm(scalarField& values) const override;

    virtual void terminalActivationTime(scalarField& values) const override;

    virtual void setTerminalActivationTime(const scalarField& values) override;

    virtual void setTerminalCoupling
    (
        const scalarField& terminalCurrent,
        const scalarField& terminalSource
    ) override;

    scalarField& membranePotential()
    {
        return Vm1D_;
    }

    scalarField& ionicCurrent()
    {
        return Iion1D_;
    }

    ionicModel& ionicModelRef()
    {
        if (!ionicModelPtr_.valid())
        {
            FatalErrorInFunction
                << "This conduction-system solver requires an ionic model, "
                << "but none was constructed."
                << exit(FatalError);
        }

        return ionicModelPtr_();
    }

    const labelList& edgeStartNodes() const
    {
        return graph_.edgeNodeA;
    }

    const labelList& edgeEndNodes() const
    {
        return graph_.edgeNodeB;
    }

    const scalarField& edgeLengths() const
    {
        return graph_.edgeLengths;
    }

    const scalarField& edgeConductances() const
    {
        return graph_.edgeConductances;
    }

    scalar chi() const
    {
        return chi_;
    }

    scalar Cm() const
    {
        return Cm_;
    }

    const conductionGraph& graph() const
    {
        return graph_;
    }

    scalarField& activationTime()
    {
        return activationTime_;
    }

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
