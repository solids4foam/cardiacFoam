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

#include "IOstreams.H"
#include "error.H"
#include "ionicSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

label ionicSelector::selectTissue(const dictionary& dict,
                                  const List<word>& supportedTissues)
{
    if (!dict.found("tissue"))
    {
        FatalErrorInFunction
            << "Normal model requires 'tissue' entry"
            << exit(FatalError);
    }

    word t;
    dict.lookup("tissue") >> t;

    if (!supportedTissues.contains(t))
    {
        FatalErrorInFunction
            << "Unsupported tissue '" << t
            << "' allowed: " << supportedTissues
            << exit(FatalError);
    }

    return tissueFlag(t);
}


//- Map tissue name to flag integer
label ionicSelector::tissueFlag(const word& t)
{
    return (t == "epicardialCells")    ? 1
         : (t == "mCells")            ? 2
         : (t == "endocardialCells")  ? 3
         : (t == "myocyte")           ? 4
                                      : -1;
}


label ionicSelector::selectDimension(const dictionary& dict,
                                     const List<word>& supportedDimensions)
{
    if (!dict.found("dimension"))
    {
        FatalErrorInFunction
            << "Model requires 'dimension' entry"
            << exit(FatalError);
    }

    word d;
    dict.lookup("dimension") >> d;

    if (!supportedDimensions.contains(d))
    {
        FatalErrorInFunction
            << "Unsupported dimension '" << d
            << "' allowed: " << supportedDimensions
            << exit(FatalError);
    }

    return dimensionFlag(d);
}


//- Map dimension name to flag integer
label ionicSelector::dimensionFlag(const word& d)
{
    return (d == "1D") ? 1 : (d == "2D") ? 2 : (d == "3D") ? 3 : -1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
