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
    Restitution curves mapping diastolic interval (DI) to action potential
    duration (APD) and conduction velocity (CV). Both curves are tabulated
    pairs read from a dictionary and evaluated by clamped piecewise-linear
    interpolation. Dimension-agnostic: shared by 1D and 3D dynamic eikonal
    solvers.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    restitutionModel.C
\*---------------------------------------------------------------------------*/

#ifndef restitutionModel_H
#define restitutionModel_H

#include "dictionary.H"
#include "scalarList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class restitutionModel Declaration
\*---------------------------------------------------------------------------*/

class restitutionModel
{
    scalarList apdDI_;
    scalarList apdVal_;
    scalarList cvDI_;
    scalarList cvVal_;

    scalar DImin_;
    scalar DImax_;
    scalar APDmin_;
    scalar APDmax_;
    scalar CVmin_;
    scalar CVmax_;

    // Read a (DI value) pair list, sorted-ascending in DI, into x and y.
    static void readCurve
    (
        const dictionary& dict,
        const word& key,
        scalarList& x,
        scalarList& y
    );

    // Clamped piecewise-linear interpolation; x must be ascending.
    static scalar interpolate
    (
        const scalarList& x,
        const scalarList& y,
        const scalar xq
    );

public:

    restitutionModel();

    //- APD for a given diastolic interval (input and output clamped).
    scalar apd(const scalar DI) const;

    //- CV for a given diastolic interval (input and output clamped).
    scalar cv(const scalar DI) const;

    //- Lower DI clamp (used to flag short-DI / railed-curve activations).
    scalar diMin() const noexcept
    {
        return DImin_;
    }

    //- Upper DI clamp (used to sanitise diagnostic DI fields).
    scalar diMax() const noexcept
    {
        return DImax_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
