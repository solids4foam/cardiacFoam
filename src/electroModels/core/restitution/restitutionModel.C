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

#include "restitutionModel.H"

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::restitutionModel::readCurve
(
    const dictionary& dict,
    const word& key,
    scalarList& x,
    scalarList& y
)
{
    const List<scalarList> entries(dict.get<List<scalarList>>(key));

    if (entries.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Restitution curve '" << key << "' is empty."
            << exit(FatalIOError);
    }

    x.setSize(entries.size());
    y.setSize(entries.size());

    forAll(entries, i)
    {
        const scalarList& e = entries[i];

        if (e.size() != 2)
        {
            FatalIOErrorInFunction(dict)
                << "Each entry in '" << key << "' must be a (DI value) pair. "
                << "Entry " << i << " has " << e.size() << " values."
                << exit(FatalIOError);
        }

        x[i] = e[0];
        y[i] = e[1];

        if (i > 0 && x[i] <= x[i - 1])
        {
            FatalIOErrorInFunction(dict)
                << "Entries in '" << key << "' must be strictly ascending in "
                << "DI. Entry " << i << " (DI=" << x[i] << ") does not follow "
                << "entry " << (i - 1) << " (DI=" << x[i - 1] << ")."
                << exit(FatalIOError);
        }
    }
}


Foam::scalar Foam::restitutionModel::interpolate
(
    const scalarList& x,
    const scalarList& y,
    const scalar xq
)
{
    if (xq <= x.first())
    {
        return y.first();
    }

    if (xq >= x.last())
    {
        return y.last();
    }

    label hi = 1;
    while (xq > x[hi])
    {
        ++hi;
    }

    const scalar w = (xq - x[hi - 1])/(x[hi] - x[hi - 1]);

    return y[hi - 1] + w*(y[hi] - y[hi - 1]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::restitutionModel::restitutionModel(const dictionary& dict)
:
    DImin_(dict.get<scalar>("DI_min")),
    DImax_(dict.get<scalar>("DI_max")),
    APDmin_(dict.get<scalar>("APD_min")),
    APDmax_(dict.get<scalar>("APD_max")),
    CVmin_(dict.get<scalar>("CV_min")),
    CVmax_(dict.get<scalar>("CV_max"))
{
    readCurve(dict, "APD_curve", apdDI_, apdVal_);
    readCurve(dict, "CV_curve", cvDI_, cvVal_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::restitutionModel::apd(const scalar DI) const
{
    const scalar DIc = min(max(DI, DImin_), DImax_);
    const scalar value = interpolate(apdDI_, apdVal_, DIc);

    return min(max(value, APDmin_), APDmax_);
}


Foam::scalar Foam::restitutionModel::cv(const scalar DI) const
{
    const scalar DIc = min(max(DI, DImin_), DImax_);
    const scalar value = interpolate(cvDI_, cvVal_, DIc);

    return min(max(value, CVmin_), CVmax_);
}


// ************************************************************************* //
