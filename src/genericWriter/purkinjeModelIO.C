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

#include "purkinjeModelIO.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

autoPtr<OFstream> purkinjeModelIO::openTimeSeries
(
    const fileName& outDir,
    const word& filename,
    const wordList& columnNames
)
{
    mkDir(outDir);
    autoPtr<OFstream> osPtr(new OFstream(outDir/filename));

    OFstream& os = osPtr.ref();
    os.setf(std::ios::scientific);
    os.precision(8);

    os << "# time";
    forAll(columnNames, i)
    {
        os << "  " << columnNames[i];
    }
    os << nl;

    return osPtr;
}


void purkinjeModelIO::writeRow
(
    OFstream& os,
    scalar time,
    const List<scalar>& values
)
{
    os << time;
    forAll(values, i)
    {
        os << "  " << values[i];
    }
    os << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
