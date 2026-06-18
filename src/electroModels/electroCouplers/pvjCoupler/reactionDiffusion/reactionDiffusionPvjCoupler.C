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


void reactionDiffusionPvjCoupler::evaluateCoupling(const char* phaseName) const
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
    R_pvj_(),
    debugCoupling_(dict.lookupOrDefault<Switch>("debugCoupling", false)),
    tissueVmBuffer_(),
    networkVmBuffer_()
{
    const scalarField* pRes = networkTerminalDomain_.terminalResistances();
    if (pRes)
    {
        R_pvj_ = *pRes;
    }
    else
    {
        R_pvj_ = scalarField(networkTerminalDomain_.terminalNodes().size(), dict.get<scalar>("rPvj"));
    }

    if (reportSetup_)
    {
        Info << "Reaction-diffusion PVJ coupling model: R_pvj=" << dict.get<scalar>("rPvj")
             << ", pvjRadius=" << pvjRadius_
             << ", couplingMode=" << couplingModeName(couplingMode_)
             << ", debugCoupling=" << debugCoupling_
             << endl;
    }
}


void reactionDiffusionPvjCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    pvjCoupler::prepareSecondaryCoupling(t0, dt);

    evaluateCoupling("secondary");

    if (couplingMode_ == unidirectional)
    {
        clearTerminalCouplingBuffers();
    }

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}


void reactionDiffusionPvjCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    evaluateCoupling("primary");

    mapper_.depositCoupling
    (
        terminalCurrentBuffer_,
        primaryDomain_.sourceField()
    );

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}

} // End namespace Foam

// ************************************************************************* //
