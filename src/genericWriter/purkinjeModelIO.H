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
    Foam::purkinjeModelIO

Description
    Static utility class for Purkinje network model I/O operations.

    Provides file-writing helpers used by purkinjeNetworkModel:

      - openTimeSeries : create a column-header .dat file
      - writeRow       : append one timestep row

    The class has no member data; all methods are static.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    purkinjeModelIO.C
\*---------------------------------------------------------------------------*/

#ifndef purkinjeModelIO_H
#define purkinjeModelIO_H

#include "OFstream.H"
#include "autoPtr.H"
#include "wordList.H"
#include "scalarList.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class purkinjeModelIO Declaration
\*---------------------------------------------------------------------------*/

class purkinjeModelIO
{
public:

    // Constructors

        //- Not instantiable – all methods are static
        purkinjeModelIO() = delete;


    // Static Member Functions

        //- Open/create a time-series output file and write the column header.
        //  Creates outDir if it does not exist.
        //  Returns an autoPtr<OFstream> to the open file.
        static autoPtr<OFstream> openTimeSeries
        (
            const fileName& outDir,
            const word& filename,
            const wordList& columnNames
        );

        //- Append one data row: time followed by all values
        static void writeRow
        (
            OFstream& os,
            scalar time,
            const List<scalar>& values
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
