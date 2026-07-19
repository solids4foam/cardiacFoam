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

#include "staggeredElectrophysicsAdvanceScheme.H"
#include "electrophysicsSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "clockTime.H"

namespace Foam
{

// Terminology used by the orchestration layer:
//   primary domain        = myocardium
//   conduction domains    = pre-myocardium graph / Purkinje domains
//   ECG domains           = post-myocardium ECG domains

defineTypeNameAndDebug(staggeredElectrophysicsAdvanceScheme, 0);
addToRunTimeSelectionTable
(
    electrophysicsAdvanceScheme,
    staggeredElectrophysicsAdvanceScheme,
    dictionary
);


staggeredElectrophysicsAdvanceScheme::
staggeredElectrophysicsAdvanceScheme(const dictionary& dict)
:
    bathPdeCouplingMethod_
    (
        dict.lookupOrDefault<word>("bathPdeCouplingMethod", "onePass")
    )
{
    if
    (
        bathPdeCouplingMethod_ != "onePass"
     && bathPdeCouplingMethod_ != "predictorCorrector"
    )
    {
        FatalErrorInFunction
            << "Unknown bathPdeCouplingMethod '"
            << bathPdeCouplingMethod_ << "'. Valid values are onePass and "
            << "predictorCorrector."
            << exit(FatalError);
    }
}

bool staggeredElectrophysicsAdvanceScheme::advance
(
    scalar t0,
    scalar dt,
    electrophysicsSystem& system,
    pimpleControl* pimplePtr,
    electrophysicsAdvanceTimings& timings
)
{
    timings = electrophysicsAdvanceTimings();

    myocardiumDomainInterface& myocardium = system.myocardium();

    clockTime timer;

    myocardium.prepareTimeStep(t0, dt);

    system.prepareConductionCouplings(t0, dt);
    timings.couplingTime = timer.timeIncrement();

    system.advanceConductionDomains(t0, dt);
    timings.conductionDomainTime = timer.timeIncrement();

    system.prepareMyocardiumCouplings(t0, dt);
    timings.couplingTime += timer.timeIncrement();

    if (system.hasPotentialDomain())
    {
        if (!myocardium.supportsSplitReactionDiffusion())
        {
            FatalErrorInFunction
                << "Unified phiE requires a split reaction/diffusion "
                << "myocardium domain."
                << exit(FatalError);
        }

        myocardium.solveReactionStep(t0, dt);
        system.preparePotentialDomain(t0, dt);

        if (bathPdeCouplingMethod_ == "onePass")
        {
            // Historical algorithm: solve phiE once from Vm at the start of
            // the step, then apply the configured PIMPLE loops only to Vm.
            system.advancePotentialDomain(t0, dt);
            myocardium.solveDiffusionStep(t0, dt, pimplePtr);
        }
        else if (bathPdeCouplingMethod_ == "predictorCorrector")
        {
            // Predict Vm with the retained phiE from the previous step,
            // update phiE from that predictor, then correct Vm. Both Vm
            // solves refer to the same old-time field and therefore do not
            // advance physical time twice.
            myocardium.solveDiffusionStepOnce(t0, dt, pimplePtr);
            system.advancePotentialDomain(t0, dt);
            myocardium.solveDiffusionStepOnce(t0, dt, pimplePtr);
        }

        myocardium.finalizeDiffusionStep();
    }
    else
    {
        myocardium.advance(t0, dt, pimplePtr);
    }
    timings.primaryDomainTime = timer.timeIncrement();

    system.prepareECGCouplings(t0, dt);
    timings.ecgCouplingTime = timer.timeIncrement();

    system.advanceECGDomains(t0, dt);
    timings.ecgDomainTime = timer.timeIncrement();

    return true;
}

} // End namespace Foam

// ************************************************************************* //
