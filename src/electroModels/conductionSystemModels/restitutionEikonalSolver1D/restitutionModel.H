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
    Foam::restitutionModel

Description
    Conduction velocity restitution curve mapping diastolic interval (DI) to
    CV. Tabulated from restitutionTemplates and evaluated by clamped
    piecewise-linear interpolation.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    restitutionModel.C
\*---------------------------------------------------------------------------*/

#ifndef restitutionModel_H
#define restitutionModel_H

#include "scalarList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class restitutionModel Declaration
\*---------------------------------------------------------------------------*/

class restitutionModel
{
    scalarList cvDI_;
    scalarList cvVal_;

    scalar DImin_;
    scalar DImax_;
    scalar CVmin_;
    scalar CVmax_;

    static scalar interpolate
    (
        const scalarList& x,
        const scalarList& y,
        const scalar xq
    );

public:

    restitutionModel();

    scalar cv(const scalar DI) const;

    scalar diMin() const noexcept
    {
        return DImin_;
    }

    scalar diMax() const noexcept
    {
        return DImax_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
