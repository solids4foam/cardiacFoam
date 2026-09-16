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

#include "eikonalMonodomainPvjCoupler.H"
#include "restitutionTemplates.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(eikonalMonodomainPvjCoupler, 0);
addToRunTimeSelectionTable
(
    electroDomainCoupler,
    eikonalMonodomainPvjCoupler,
    dictionary
);


eikonalMonodomainPvjCoupler::eikonalMonodomainPvjCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    pvjCoupler(primaryDomain, secondaryDomain, dict),
    terminalActivationBuffer_(),
    R_pvj_(),
    vmTemplateOffset_(0.0)
{
    if (dict.found("rPvj"))
    {
        R_pvj_ = scalarField(networkTerminalDomain_.terminalNodes().size(), dict.get<scalar>("rPvj"));
    }
    else
    {
        FatalErrorInFunction
            << "Missing rPvj in domainCouplings dictionary for "
            << "eikonalMonodomainPvjCoupler" << exit(FatalError);
    }

    const ionicModel* tissueModel = primaryDomain.ionicModelPtr();
    if (!tissueModel)
    {
        FatalErrorInFunction
            << "eikonalMonodomainPvjCoupler requires a tissue domain with an "
            << "ionic model" << exit(FatalError);
    }

    // Resting potential from ranks that hold tissue cells.
    const PtrList<scalarField>* tissueStates = tissueModel->ioStatesPtr();
    const bool hasCells =
        tissueStates && !tissueStates->empty() && !(*tissueStates)[0].empty();
    scalar vmRestSum = hasCells ? tissueModel->vmRest() : 0.0;
    label nRanksWithCells = hasCells ? 1 : 0;
    reduce(vmRestSum, sumOp<scalar>());
    reduce(nRanksWithCells, sumOp<label>());

    vmTemplateOffset_ =
        vmRestSum/max(nRanksWithCells, 1)
      - restitutionTemplates::purkinjeVmValues[0]*1e-3;
}


void eikonalMonodomainPvjCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    if (couplingMode_ == bidirectional)
    {
        observeTerminalActivations();
    }

    clearTerminalCouplingBuffers();

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}


void eikonalMonodomainPvjCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    networkTerminalDomain_.terminalActivationTime(terminalActivationBuffer_);

    // Anterograde coupling using voltage template
    const label nTerminalNodes = networkTerminalDomain_.terminalNodes().size();

    // Template shifted to the tissue resting potential, so a quiescent node
    // drives no junction current and leaves the junction tissue at rest.
    scalarField terminalVoltage
    (
        nTerminalNodes,
        restitutionTemplates::purkinjeVmValues[0]*1e-3 + vmTemplateOffset_
    );
    const scalar currentTime = mesh_.time().value();

    for (label i = 0; i < nTerminalNodes; ++i)
    {
        const scalar tact = terminalActivationBuffer_[i];

        // If activated, compute the template voltage
        if (currentTime >= tact)
        {
            const scalar localTime = currentTime - tact;
            terminalVoltage[i] =
                restitutionTemplates::evaluatePurkinjeVmTemplate(localTime)*1e-3
              + vmTemplateOffset_;
        }
    }

    scalarField tissueVoltage;
    mapper_.gatherVm3DPvjs(primaryDomain_.Vm(), tissueVoltage);

    const scalarField* specificR_pvj = networkTerminalDomain_.terminalResistances();
    if (specificR_pvj)
    {
        R_pvj_ = *specificR_pvj;
    }

    for (label i = 0; i < nTerminalNodes; ++i)
    {
        terminalCurrentBuffer_[i] = (terminalVoltage[i] - tissueVoltage[i]) / R_pvj_[i];
    }

    mapper_.depositCoupling
    (
        terminalCurrentBuffer_,
        primaryDomain_.sourceField()
    );

    mapper_.volumetricSource
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}

} // End namespace Foam

// ************************************************************************* //
