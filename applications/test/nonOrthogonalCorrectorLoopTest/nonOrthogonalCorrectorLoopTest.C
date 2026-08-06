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

#include "nonOrthogonalCorrectorLoop.H"
#include "IOstreams.H"
#include <cassert>

using namespace Foam;

int main()
{
    // count overload: N correctors => N+1 assemble-and-solve passes
    for (const label n : {label(0), label(1), label(2), label(5)})
    {
        label calls = 0;
        correctNonOrthogonalLoop(n, [&]() { ++calls; });
        assert(calls == n + 1);
    }

    // negative count clamps to a single pass
    {
        label calls = 0;
        correctNonOrthogonalLoop(label(-3), [&]() { ++calls; });
        assert(calls == 1);
    }

    Info<< "nonOrthogonalCorrectorLoop count semantics OK" << endl;
    return 0;
}
