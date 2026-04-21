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

#ifndef solverStartupLog_H
#define solverStartupLog_H

#include "IOstreams.H"
#include "word.H"

namespace Foam
{
namespace solverStartupLog
{

inline void readingDictionary(const word& dictName)
{
    Info<< "Reading " << dictName << nl << endl;
}

inline void sectionDivider()
{
    Info<< "---------------------------------------------" << nl;
}

inline void logSolverMode(const bool useExplicitAlgorithm)
{
    Info<< nl
        << "Using " << (useExplicitAlgorithm ? "explicit" : "implicit")
        << " solver" << nl << endl;
}

} // End namespace solverStartupLog
} // End namespace Foam

#endif
