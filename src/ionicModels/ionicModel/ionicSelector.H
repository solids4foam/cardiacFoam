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

#ifndef ionicSelector_H
#define ionicSelector_H

#include "List.H"
#include "dictionary.H"
#include "word.H"

namespace Foam {

class ionicSelector {
public:
  // Main interface
  static label selectTissue(const dictionary &dict,
                            const List<word> &supportedTissues);
  static label selectDimension(const dictionary &dict,
                               const List<word> &supportedDimensions);

  static label tissueFlag(const word &t);
  static label dimensionFlag(const word &d);
};

} // End namespace Foam

#endif
