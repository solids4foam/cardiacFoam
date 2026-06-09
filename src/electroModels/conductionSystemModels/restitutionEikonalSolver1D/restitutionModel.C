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

#include "restitutionTemplates.H"

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

Foam::restitutionModel::restitutionModel()
:
    DImin_(restitutionTemplates::purkinjeDI[0]),
    DImax_(restitutionTemplates::purkinjeDI[restitutionTemplates::purkinjeSize - 1]),
    APDmin_(restitutionTemplates::purkinjeAPD[0]),
    APDmax_(restitutionTemplates::purkinjeAPD[restitutionTemplates::purkinjeSize - 1]),
    CVmin_(restitutionTemplates::purkinjeCV[0]),
    CVmax_(restitutionTemplates::purkinjeCV[restitutionTemplates::purkinjeSize - 1])
{
    using namespace restitutionTemplates;

    apdDI_.setSize(purkinjeSize);
    apdVal_.setSize(purkinjeSize);
    for(label i = 0; i < purkinjeSize; ++i)
    {
        apdDI_[i] = purkinjeDI[i];
        apdVal_[i] = purkinjeAPD[i];
    }

    cvDI_.setSize(purkinjeSize);
    cvVal_.setSize(purkinjeSize);
    for(label i = 0; i < purkinjeSize; ++i)
    {
        cvDI_[i] = purkinjeDI[i];
        cvVal_[i] = purkinjeCV[i];
    }
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
