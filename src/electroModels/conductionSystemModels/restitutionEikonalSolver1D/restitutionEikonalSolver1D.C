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

#include "restitutionEikonalSolver1D.H"
#include "conductionSystemDomain.H"
#include "restitutionTemplates.H"
#include "addToRunTimeSelectionTable.H"
#include <queue>
#include <utility>
#include <vector>

namespace Foam
{

defineTypeNameAndDebug(restitutionEikonalSolver1D, 0);
addToRunTimeSelectionTable
(
    conductionSystemSolver,
    restitutionEikonalSolver1D,
    dictionary
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::restitutionEikonalSolver1D::restitutionEikonalSolver1D
(
    const fvMesh&,
    const dictionary& solverCoeffs
)
:
    restitutionPtr_(new restitutionModel()),
    apdNominal_
    (
        solverCoeffs.lookupOrDefault<scalar>
        (
            "apdNominal",
            restitutionTemplates::purkinjeAPDnominal
        )
    ),
    minBeatInterval_(apdNominal_ + restitutionPtr_->diMin()),
    escapeInterval_
    (
        solverCoeffs.lookupOrDefault<scalar>
        (
            "escapeInterval",
            restitutionTemplates::purkinjeEscapeInterval
        )
    ),
    tStart_(0),
    useEdgeConductance_
    (
        solverCoeffs.lookupOrDefault<Switch>("useEdgeConductance", true)
    ),
    refConductance_
    (
        solverCoeffs.lookupOrDefault<scalar>("referenceConductance", 1.0)
    ),
    initialised_(false)
{
    if (solverCoeffs.found("stimulus"))
    {
        const dictionary& sDict = solverCoeffs.subDict("stimulus");

        stimSites_ = sDict.get<labelList>("sites");
        stimProtocol_ = stimulusIO::loadStimulusProtocol(sDict);
    }
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::restitutionEikonalSolver1D::initialiseState
(
    conductionSystemDomain& domain,
    const scalar t0
)
{
    const label N = domain.graph().nNodes;
    const scalarField& Tact = domain.activationTime();

    lastActTime_.setSize(N, -GREAT);
    DI_.setSize(N, GREAT);
    nextTact_.setSize(N, GREAT);
    minDI_.setSize(N, GREAT);

    activatedFrom_.setSize(N, -1);
    nextTactSource_.setSize(N, -1);

    blockCount_.setSize(N, 0);
    wavebreakCount_.setSize(N, 0);

    tStart_ = t0;

    // Seed the event queue from any activation times already present on the
    // graph (for example the rootNode activation configured at t=0, or a
    // restart that carries past/future graph activations).
    forAll(Tact, i)
    {
        if (Tact[i] < 0.0)
        {
            continue;
        }

        if (Tact[i] >= t0 - SMALL)
        {
            nextTact_[i] = Tact[i];
            nextTactSource_[i] = -1;
        }
        else
        {
            lastActTime_[i] = Tact[i];
        }
    }

    initialised_ = true;
}


void Foam::restitutionEikonalSolver1D::importExternalActivations
(
    conductionSystemDomain& domain,
    const scalar tNow
)
{
    scalarField& Tact = domain.activationTime();
    const labelList& terminalNodes = domain.terminalNodes();

    forAll(terminalNodes, i)
    {
        const label nodeI = terminalNodes[i];
        const scalar incomingTime = Tact[nodeI];
        const scalar lastActTime = lastActTime_[nodeI];

        // Negative activation times mean "not yet activated" for this graph
        // field and must not be treated as a real incoming wave.
        if (incomingTime < 0.0)
        {
            continue;
        }

        if (incomingTime <= lastActTime + SMALL)
        {
            continue;
        }

        Tact[nodeI] = lastActTime < 0.0 ? -1.0 : lastActTime;

        if (incomingTime > tNow)
        {
            if (incomingTime < nextTact_[nodeI])
            {
                nextTact_[nodeI] = incomingTime;
                nextTactSource_[nodeI] = -1;
            }

            continue;
        }

        if
        (
            incomingTime - lastActTime >= minBeatInterval_
        )
        {
            if (incomingTime < nextTact_[nodeI])
            {
                nextTact_[nodeI] = incomingTime;
                nextTactSource_[nodeI] = -1;
            }
        }
        else
        {
            ++blockCount_[nodeI];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::restitutionEikonalSolver1D::advance
(
    conductionSystemDomain& domain,
    scalar t0,
    scalar dt
)
{
    if (!initialised_)
    {
        initialiseState(domain, t0);
    }

    const scalar tNow = t0 + dt;

    const conductionGraph& G = domain.graph();
    scalarField& Tact = domain.activationTime();

    importExternalActivations(domain, tNow);

    forAll(lastActTime_, i)
    {
        const scalar tRef =
            lastActTime_[i] < 0 ? tStart_ : lastActTime_[i];
        const scalar tEscape = tRef + escapeInterval_;
        if (tEscape < nextTact_[i])
        {
            nextTact_[i] = tEscape;
            nextTactSource_[i] = -1;
        }
    }

    if (stimulusIO::computeStimulus(tNow, stimProtocol_) != 0)
    {
        forAll(stimSites_, s)
        {
            const label site = stimSites_[s];
            const scalar beatInterval =
                lastActTime_[site] < 0 ? GREAT : tNow - lastActTime_[site];

            if (beatInterval >= minBeatInterval_ && tNow < nextTact_[site])
            {
                nextTact_[site] = tNow;
                nextTactSource_[site] = -1;
            }
        }
    }

    using Event = std::pair<scalar, label>;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> pq;

    forAll(nextTact_, i)
    {
        if (nextTact_[i] <= tNow)
        {
            const scalar beatInterval =
                lastActTime_[i] < 0 ? GREAT : nextTact_[i] - lastActTime_[i];
            if (beatInterval >= minBeatInterval_)
            {
                pq.push(std::make_pair(nextTact_[i], i));
            }
        }
    }

    while (!pq.empty())
    {
        const scalar te = pq.top().first;
        const label  i  = pq.top().second;
        pq.pop();

        const scalar beatInterval_i =
            lastActTime_[i] < 0 ? GREAT : te - lastActTime_[i];
        if (beatInterval_i < minBeatInterval_ || te != nextTact_[i])
        {
            continue;
        }

        const scalar DIact = beatInterval_i - apdNominal_;

        Tact[i] = te;
        DI_[i] = DIact;
        lastActTime_[i] = te;

        activatedFrom_[i] = nextTactSource_[i];
        nextTact_[i] = GREAT;
        nextTactSource_[i] = -1;

        if (DIact < minDI_[i])
        {
            minDI_[i] = DIact;
        }

        const label from = activatedFrom_[i];

        for (label k = G.adjOffsets[i]; k < G.adjOffsets[i + 1]; ++k)
        {
            const label j = G.adjNeighbours[k];

            if (j == from)
            {
                continue;
            }

            const scalar beatInterval_j =
                lastActTime_[j] < 0 ? GREAT : te - lastActTime_[j];
            if (beatInterval_j < minBeatInterval_)
            {
                ++blockCount_[j];
                ++wavebreakCount_[i];
                continue;
            }

            const label eI = G.adjEdges[k];
            const scalar DIj = beatInterval_j - apdNominal_;
            scalar cv = restitutionPtr_->cv(DIj);

            if (useEdgeConductance_)
            {
                const scalar gRel = G.edgeConductances[eI]/refConductance_;

                if (gRel <= SMALL)
                {
                    ++blockCount_[j];
                    ++wavebreakCount_[i];
                    continue;
                }

                cv *= sqrt(gRel);
            }

            const scalar cand = te + G.edgeLengths[eI]/cv;

            if (cand < nextTact_[j])
            {
                nextTact_[j] = cand;
                nextTactSource_[j] = i;

                if (cand <= tNow)
                {
                    pq.push(std::make_pair(cand, j));
                }
            }
        }
    }

    scalarField& Vm = domain.membranePotential();
    forAll(Vm, i)
    {
        if (Tact[i] < 0.0)
        {
            Vm[i] = restitutionTemplates::purkinjeVmValues[0] * 1e-3;
        }
        else
        {
            const scalar localTime = tNow - Tact[i];
            Vm[i] = restitutionTemplates::evaluatePurkinjeVmTemplate(localTime) * 1e-3;
        }
    }
}


void Foam::restitutionEikonalSolver1D::diagnosticFields
(
    wordList& names,
    PtrList<scalarField>& fields
) const
{
    const label N = lastActTime_.size();

    if (N == 0)
    {
        return;
    }

    const scalar diMaxV = restitutionPtr_->diMax();

    names.setSize(4);
    fields.setSize(4);

    scalarField* diPtr = new scalarField(N, 0.0);
    forAll(*diPtr, i)
    {
        (*diPtr)[i] = min(max(DI_[i], scalar(0)), diMaxV);
    }
    names[0] = "DI";
    fields.set(0, diPtr);

    scalarField* minDIPtr = new scalarField(N, 0.0);
    forAll(*minDIPtr, i)
    {
        (*minDIPtr)[i] =
        (
            minDI_[i] >= GREAT
          ? diMaxV
          : min(max(minDI_[i], scalar(0)), diMaxV)
        );
    }
    names[1] = "minDI";
    fields.set(1, minDIPtr);

    scalarField* blockPtr = new scalarField(N, 0.0);
    forAll(*blockPtr, i)
    {
        (*blockPtr)[i] = scalar(blockCount_[i]);
    }
    names[2] = "blockCount";
    fields.set(2, blockPtr);

    scalarField* wbPtr = new scalarField(N, 0.0);
    forAll(*wbPtr, i)
    {
        (*wbPtr)[i] = scalar(wavebreakCount_[i]);
    }
    names[3] = "wavebreakCount";
    fields.set(3, wbPtr);
}


// ************************************************************************* //
