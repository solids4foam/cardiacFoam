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

#include "reactionDiffusionPvjCoupler.H"

#include "electroVolumeFieldDomain.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(reactionDiffusionPvjCoupler, 0);
addToRunTimeSelectionTable
(
    electroDomainCoupler,
    reactionDiffusionPvjCoupler,
    dictionary
);


void reactionDiffusionPvjCoupler::couplingCurrentAtPvjs
(
    const scalarField& networkVm,
    const scalarField& tissueVm,
    scalarField& current
) const
{
    if (networkVm.size() != tissueVm.size())
    {
        FatalErrorInFunction
            << "PVJ coupling expected matching field sizes but received "
            << networkVm.size() << " and " << tissueVm.size()
            << exit(FatalError);
    }

    current.setSize(networkVm.size());
    forAll(current, i)
    {
        current[i] = (networkVm[i] - tissueVm[i])/R_pvj_[i];
    }
}


void reactionDiffusionPvjCoupler::evaluateCoupling
(
    scalar primaryTime,
    scalar secondaryTime,
    const char* phaseName
)
{
    mapper_.gatherVm3DPvjs(primaryDomain_.Vm(), tissueVmBuffer_);
    networkTerminalDomain_.terminalVm(networkVmBuffer_);
    couplingCurrentAtPvjs
    (
        networkVmBuffer_,
        tissueVmBuffer_,
        terminalCurrentBuffer_
    );
    mapper_.volumetricSource(terminalCurrentBuffer_, terminalSourceBuffer_);

    reportCouplingDiagnostics(phaseName);
}


void reactionDiffusionPvjCoupler::reportCouplingDiagnostics
(
    const char* phaseName
) const
{
    if (!debugCoupling_ || terminalCurrentBuffer_.empty())
    {
        return;
    }

    Info<< "PVJ coupling debug (" << phaseName << "): "
        << "networkVm[min,max]=[" << gMin(networkVmBuffer_)
        << ", " << gMax(networkVmBuffer_) << "], "
        << "tissueVm[min,max]=[" << gMin(tissueVmBuffer_)
        << ", " << gMax(tissueVmBuffer_) << "], "
        << "terminalCurrent[min,max]=[" << gMin(terminalCurrentBuffer_)
        << ", " << gMax(terminalCurrentBuffer_) << "], "
        << "terminalSource[min,max]=[" << gMin(terminalSourceBuffer_)
        << ", " << gMax(terminalSourceBuffer_) << "]"
        << nl;
}


reactionDiffusionPvjCoupler::reactionDiffusionPvjCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    pvjCoupler(primaryDomain, secondaryDomain, dict),
    R_pvj_(readResistances(dict)),
    debugCoupling_(dict.lookupOrDefault<Switch>("debugCoupling", false)),
    couplingScheme_(readCouplingScheme(dict)),
    tissueVmBuffer_(),
    networkVmBuffer_()
{
    if (couplingMode_ == bidirectional)
    {
        networkTerminalDomain_.setTerminalConductances(1.0/R_pvj_);

        if (couplingScheme_ == "implicit")
        {
            warnImplicitChargeLag();
        }
    }
}


void reactionDiffusionPvjCoupler::warnImplicitChargeLag() const
{
    // The network loses G*(Vn' - <V>) against the tissue average before the
    // tissue solve, the tissue gains G*(Vn' - <V>') after it, so each step
    // misplaces G*dt*(<V>' - <V>): a spurious capacitance dt/R per junction.
    const electroVolumeFieldDomain* tissue =
        dynamic_cast<const electroVolumeFieldDomain*>(&primaryDomain_);

    if (!tissue)
    {
        return;
    }

    const scalar dt = mesh_.time().deltaTValue();
    const scalar chiCm = (tissue->chi()*tissue->Cm()).value();
    const scalarField ratio
    (
        dt/(R_pvj_*chiCm*mapper_.sphereVolumes())
    );

    WarningInFunction
        << "pvjCouplingScheme implicit with couplingMode bidirectional "
        << "does not conserve charge per step: the staggered solve acts as "
        << "an extra capacitance deltaT/R at each junction, up to "
        << gMax(ratio) << " times the tissue's chi*Cm*V_s at deltaT = "
        << dt << " s. pvjCouplingScheme explicit conserves it exactly."
        << endl;
}


void reactionDiffusionPvjCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    pvjCoupler::prepareSecondaryCoupling(t0, dt);

    evaluateCoupling(t0, t0, "secondary");

    // The network solves its junction term at t0 + dt against the tissue
    // at t0; the manufactured source matches those levels.
    if (verificationModelPtr_)
    {
        verificationModelPtr_->updateManufacturedSource
        (
            primaryDomain_,
            secondaryDomain_,
            t0,
            t0 + dt,
            couplingScheme_ == "implicit",
            couplingMode_ == bidirectional,
            "secondary"
        );
    }

    if (couplingMode_ == bidirectional)
    {
        networkTerminalDomain_.setTerminalTissueVm(tissueVmBuffer_);
    }
}


void reactionDiffusionPvjCoupler::depositPrimaryCoupling(const scalar dt) const
{
    if (couplingScheme_ == "explicit")
    {
        checkExplicitCouplingStability(R_pvj_, dt);

        mapper_.depositCoupling
        (
            terminalCurrentBuffer_,
            primaryDomain_.sourceField()
        );
        return;
    }

    volScalarField* implicitSourceCoeff =
        primaryDomain_.implicitSourceCoeffPtr();

    if (!implicitSourceCoeff)
    {
        FatalErrorInFunction
            << "pvjCouplingScheme implicit requires the primary tissue "
            << "domain to expose an implicit source coefficient field."
            << exit(FatalError);
    }

    mapper_.depositImplicitCoupling
    (
        networkVmBuffer_,
        R_pvj_,
        primaryDomain_.sourceField(),
        *implicitSourceCoeff
    );
}


void reactionDiffusionPvjCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    evaluateCoupling(t0, t0 + dt, "primary");

    if (verificationModelPtr_)
    {
        verificationModelPtr_->updateManufacturedSource
        (
            primaryDomain_,
            secondaryDomain_,
            t0,
            t0 + dt,
            couplingScheme_ == "implicit",
            couplingMode_ == bidirectional,
            "primary"
        );
    }

    depositPrimaryCoupling(dt);

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}

} // End namespace Foam

// ************************************************************************* //
