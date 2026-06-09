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
    useEdgeConductance_
    (
        solverCoeffs.lookupOrDefault<Switch>("useEdgeConductance", true)
    ),
    refConductance_
    (
        solverCoeffs.lookupOrDefault<scalar>("referenceConductance", 1.0)
    ),
    reportSetup_(solverCoeffs.lookupOrDefault<Switch>("reportSetup", false)),
    initialised_(false)
{
    if (solverCoeffs.found("stimulus"))
    {
        const dictionary& sDict = solverCoeffs.subDict("stimulus");

        stimSites_ = sDict.get<labelList>("sites");
        stimProtocol_ = stimulusIO::loadStimulusProtocol(sDict);
    }

    if (reportSetup_)
    {
        Info<< "restitutionEikonalSolver1D: " << stimSites_.size()
            << " stimulus sites; S1 start=" << stimProtocol_.stimStart
            << " period=" << stimProtocol_.stimPeriodS1
            << " n=" << stimProtocol_.nStim1
            << "; S2 coupling=" << stimProtocol_.stimPeriodS2
            << " n=" << stimProtocol_.nStim2 << endl;
    }
}


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::restitutionEikonalSolver1D::initialiseState
(
    const conductionSystemDomain& domain
)
{
    const label N = domain.graph().nNodes;

    RT_.setSize(N, -GREAT);
    DI_.setSize(N, GREAT);
    APD_.setSize(N, restitutionPtr_->apd(GREAT));
    nextTact_.setSize(N, GREAT);
    minDI_.setSize(N, GREAT);

    state_.setSize(N, excitable);
    activatedFrom_.setSize(N, -1);
    nextTactSource_.setSize(N, -1);

    blockCount_.setSize(N, 0);
    wavebreakCount_.setSize(N, 0);
    shortDICount_.setSize(N, 0);

    history_.reset(N);

    initialised_ = true;
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
        initialiseState(domain);
    }

    const scalar tNow = t0 + dt;

    const conductionGraph& G = domain.graph();
    scalarField& Tact = domain.activationTime();

    forAll(state_, i)
    {
        if (state_[i] == refractory && tNow >= RT_[i])
        {
            state_[i] = excitable;
        }
    }

    if (stimulusIO::computeStimulus(tNow, stimProtocol_) != 0)
    {
        forAll(stimSites_, s)
        {
            const label site = stimSites_[s];

            if (state_[site] == excitable && tNow < nextTact_[site])
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
        if (state_[i] == excitable && nextTact_[i] <= tNow)
        {
            pq.push(std::make_pair(nextTact_[i], i));
        }
    }

    while (!pq.empty())
    {
        const scalar te = pq.top().first;
        const label  i  = pq.top().second;
        pq.pop();

        if (state_[i] != excitable || te != nextTact_[i])
        {
            continue;
        }

        const scalar DIact = te - RT_[i];

        Tact[i] = te;
        history_.record(i, te);
        DI_[i] = DIact;
        APD_[i] = restitutionPtr_->apd(DIact);
        RT_[i] = te + APD_[i];

        state_[i] = refractory;
        activatedFrom_[i] = nextTactSource_[i];
        nextTact_[i] = GREAT;
        nextTactSource_[i] = -1;

        if (DIact < minDI_[i])
        {
            minDI_[i] = DIact;
        }

        if (DIact <= restitutionPtr_->diMin())
        {
            ++shortDICount_[i];
        }

        const label from = activatedFrom_[i];

        for (label k = G.adjOffsets[i]; k < G.adjOffsets[i + 1]; ++k)
        {
            const label j = G.adjNeighbours[k];

            if (j == from)
            {
                continue;
            }

            if (state_[j] == refractory)
            {
                ++blockCount_[j];
                ++wavebreakCount_[i];
                continue;
            }

            const label eI = G.adjEdges[k];
            const scalar DIj = te - RT_[j];
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
}


void Foam::restitutionEikonalSolver1D::diagnosticFields
(
    wordList& names,
    PtrList<scalarField>& fields
) const
{
    const label N = state_.size();

    if (N == 0)
    {
        return;
    }

    const scalar diMaxV = restitutionPtr_->diMax();

    names.setSize(7);
    fields.setSize(7);

    names[0] = "APD";
    fields.set(0, new scalarField(APD_));

    scalarField* diPtr = new scalarField(N, 0.0);
    forAll(*diPtr, i)
    {
        (*diPtr)[i] = min(max(DI_[i], scalar(0)), diMaxV);
    }
    names[1] = "DI";
    fields.set(1, diPtr);

    names[2] = "state";
    scalarField* statePtr = new scalarField(N, 0.0);
    forAll(*statePtr, i)
    {
        (*statePtr)[i] = scalar(state_[i]);
    }
    fields.set(2, statePtr);

    names[3] = "minDI";
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
    fields.set(3, minDIPtr);

    names[4] = "blockCount";
    scalarField* blockPtr = new scalarField(N, 0.0);
    forAll(*blockPtr, i)
    {
        (*blockPtr)[i] = scalar(blockCount_[i]);
    }
    fields.set(4, blockPtr);

    names[5] = "wavebreakCount";
    scalarField* wbPtr = new scalarField(N, 0.0);
    forAll(*wbPtr, i)
    {
        (*wbPtr)[i] = scalar(wavebreakCount_[i]);
    }
    fields.set(5, wbPtr);

    names[6] = "shortDICount";
    scalarField* sdPtr = new scalarField(N, 0.0);
    forAll(*sdPtr, i)
    {
        (*sdPtr)[i] = scalar(shortDICount_[i]);
    }
    fields.set(6, sdPtr);

    history_.appendDiagnostics(names, fields);
}


// ************************************************************************* //
