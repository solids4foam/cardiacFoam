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
    myocardiumDomainInterface

Description
    Interface for myocardium domains.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef myocardiumDomainInterface_H
#define myocardiumDomainInterface_H

#include "electroDomainInterface.H"
#include "electroStateProvider.H"
#include "electromechanicalSignalProvider.H"
#include "electroDomainCouplingEndpoints.H"
#include "autoPtr.H"
#include "pimpleControl.H"
#include "myocardiumMeshSubset.H"

namespace Foam
{

class fvMesh;
class dictionary;
class ionicModel;
class electroVerificationModel;

class myocardiumDomainInterface
:
    public electroDomainInterface,
    public tissueCouplingEndpoint,
    public electroStateProvider
{
public:

    using electroDomainInterface::advance;

    //- Build the domain of the myocardiumSolver solverType:
    //  eikonalMyocardiumDomain for its solver name, else a reaction-diffusion
    //  myocardiumDomain with that myocardiumSolver kernel. externalPhiE: the
    //  bath domain solves the extracellular potential (torsoCoupling bath).
    static autoPtr<myocardiumDomainInterface> New
    (
        const fvMesh& mesh,
        const word& solverType,
        const dictionary& electroProperties,
        scalar initialDeltaT,
        const bool externalPhiE
    );

    virtual ~myocardiumDomainInterface() = default;

    //- Read the tissue's ionic restart state; false when there is none or
    //  the tissue has no ionic model
    virtual bool readRestartState()
    {
        return false;
    }

    //- Write the tissue's ionic restart state (none without an ionic model)
    virtual void writeRestartState() const
    {}

    //- Does the tissue model solve an extracellular potential (bidomain)?
    virtual bool hasExtracellularPotential() const
    {
        return false;
    }

    virtual const ElectromechanicalSignalProvider* provider() const
    {
        return nullptr;
    }

    virtual const volScalarField& activationTime() const = 0;

    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) = 0;

    virtual bool supportsSplitReactionDiffusion() const
    {
        return false;
    }

    virtual void solveReactionStep(scalar t0, scalar dt)
    {
        (void)t0;
        (void)dt;

        FatalErrorInFunction
            << "This myocardium domain does not support split "
            << "reaction/diffusion advance."
            << exit(FatalError);
    }

    virtual void solveDiffusionStep
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    )
    {
        (void)t0;
        (void)dt;
        (void)pimplePtr;

        FatalErrorInFunction
            << "This myocardium domain does not support split "
            << "reaction/diffusion advance."
            << exit(FatalError);
    }

    //- Solve the implicit diffusion equation once for the opt-in coupled
    //  bath-PDE predictor/corrector; it must not advance the ionic model or
    //  physical time.
    virtual void solveDiffusionStepOnce
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    )
    {
        (void)t0;
        (void)dt;
        (void)pimplePtr;

        FatalErrorInFunction
            << "This myocardium domain does not support a single diffusion "
            << "sub-step."
            << exit(FatalError);
    }

    //- Complete the split step over [t0, t0+dt] once Vm(t0+dt) is solved
    virtual void finalizeDiffusionStep(scalar t0, scalar dt)
    {
        (void)t0;
        (void)dt;
    }

    //- Read the bath domain's heart phiE (heart mesh) in the myocardium
    //  solver in place of an own phiE
    virtual void bindExternalPhiE(const volScalarField& heartPhiE)
    {
        (void)heartPhiE;

        FatalErrorInFunction
            << "The selected myocardium domain cannot bind an external phiE "
            << "field. Unified phiE requires a bidomain myocardium domain."
            << exit(FatalError);
    }

    //- Verification only: grant the bath's global phiE (base mesh) to this
    //  domain's verification model, which writes the manufactured initial
    //  values into it; heartCellMap maps heart to base cells. Nothing to do
    //  without a verification model.
    virtual void bindVerificationPhiE
    (
        volScalarField& globalPhiE,
        const labelUList& heartCellMap
    )
    {
        (void)globalPhiE;
        (void)heartCellMap;
    }

    virtual void unbindExternalPhiE()
    {
        FatalErrorInFunction
            << "The selected myocardium domain cannot unbind an external phiE "
            << "field. Unified phiE requires a bidomain myocardium domain."
            << exit(FatalError);
    }

    //- Verification model of this domain, or nullptr.
    virtual const electroVerificationModel* verificationModel() const
    {
        return nullptr;
    }

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const = 0;

    virtual bool shouldPostProcess() const = 0;

    virtual void exportStates() = 0;

    virtual void postProcess() = 0;

    virtual bool applyModelTimeControls(Time& runTime)
    {
        (void)runTime;
        return false;
    }
};


} // End namespace Foam

#endif

// ************************************************************************* //
