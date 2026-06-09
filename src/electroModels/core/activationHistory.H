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
    Foam::activationHistory

Description
    Per-node multiple-activation (multi-beat) history store. Keeps an
    activation count per node and a dynamically grown list of per-beat
    activation-time fields (activationTime_1, activationTime_2, ...), each
    initialised to a safe -1 placeholder for nodes that have not yet reached
    that beat.

    This is storage only — detection stays in the owner. Eikonal solvers call
    record() natively when an activation is accepted; voltage solvers
    (monodomain / bidomain) would call it from their Vm-threshold state
    machine. Index-space agnostic: works for 1-D graph nodes or 3-D cells.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    activationHistory.C
\*---------------------------------------------------------------------------*/

#ifndef activationHistory_H
#define activationHistory_H

#include "scalarField.H"
#include "labelList.H"
#include "PtrList.H"
#include "wordList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class activationHistory Declaration
\*---------------------------------------------------------------------------*/

class activationHistory
{
    labelList            count_;   //- activations per node
    PtrList<scalarField> times_;   //- times_[b][i]: (b+1)-th activation time
    label                maxCount_;

public:

    activationHistory()
    :
        maxCount_(0)
    {}

    //- Reset to N nodes with zero activations.
    void reset(const label N)
    {
        count_.setSize(N, 0);
        times_.clear();
        maxCount_ = 0;
    }

    //- Record an activation of node i at time t. Grows the history on demand
    //  and initialises new beat fields to -1 for all other nodes.
    void record(const label i, const scalar t)
    {
        ++count_[i];

        if (count_[i] > maxCount_)
        {
            maxCount_ = count_[i];
            times_.append(new scalarField(count_.size(), -1.0));
        }

        times_[count_[i] - 1][i] = t;
    }

    //- Highest activation count reached by any node.
    label maxCount() const noexcept
    {
        return maxCount_;
    }

    //- Append "activationCount" and "activationTime_1..N" to a diagnostic
    //  field list (the lists are extended, existing entries preserved).
    void appendDiagnostics
    (
        wordList& names,
        PtrList<scalarField>& fields
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
