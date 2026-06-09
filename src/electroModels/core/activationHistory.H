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
    initialised to -1 for nodes that have not reached that beat.

    Storage only; detection stays in the owner. Eikonal solvers call record()
    when an activation is accepted; voltage solvers (monodomain / bidomain)
    call it from their Vm-threshold detection. Index-space agnostic: works for
    1-D graph nodes or 3-D cells.

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
    labelList            count_;
    PtrList<scalarField> times_;
    label                maxCount_;

public:

    activationHistory()
    :
        maxCount_(0)
    {}

    void reset(const label N)
    {
        count_.setSize(N, 0);
        times_.clear();
        maxCount_ = 0;
    }

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

    label maxCount() const noexcept
    {
        return maxCount_;
    }

    void appendDiagnostics
    (
        wordList& names,
        PtrList<scalarField>& fields
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
