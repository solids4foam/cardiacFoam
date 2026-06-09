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

#include "activationHistory.H"
#include "word.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activationHistory::appendDiagnostics
(
    wordList& names,
    PtrList<scalarField>& fields
) const
{
    const label base = names.size();
    const label add = 1 + times_.size();

    names.setSize(base + add);
    fields.setSize(base + add);

    names[base] = "activationCount";
    scalarField* countPtr = new scalarField(count_.size(), 0.0);
    forAll(count_, i)
    {
        (*countPtr)[i] = scalar(count_[i]);
    }
    fields.set(base, countPtr);

    forAll(times_, b)
    {
        names[base + 1 + b] = "activationTime_" + name(b + 1);
        fields.set(base + 1 + b, new scalarField(times_[b]));
    }
}


// ************************************************************************* //
