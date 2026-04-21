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
    Foam::ecgModelIO

Description
    Static utility class for ECG model I/O operations.

    Provides file-writing helpers used by ecgModel subclasses:

      - openTimeSeries : create a column-header .dat file
      - writeRow       : append one timestep row
      - loadSurface    : load a torso STL and prepare VTK output directory
      - writeSurfaceVtk: write face-centred potential as VTK ASCII POLYDATA

    The class has no member data; all methods are static.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ecgModelIO.C
\*---------------------------------------------------------------------------*/

#ifndef ecgModelIO_H
#define ecgModelIO_H

#include "OFstream.H"
#include "autoPtr.H"
#include "wordList.H"
#include "pointField.H"
#include "scalarList.H"
#include "fileName.H"
#include "triSurface.H"
#include "Time.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ecgModelIO Declaration
\*---------------------------------------------------------------------------*/

class ecgModelIO
{
public:

    // Constructors

        //- Not instantiable – all methods are static
        ecgModelIO() = delete;


    // Static Member Functions

        //- Open/create a time-series output file and write the column header.
        //  Creates outDir if it does not exist.
        //  Returns an autoPtr<OFstream> to the open file.
        static autoPtr<OFstream> openTimeSeries
        (
            const fileName& outDir,
            const word& filename,
            const wordList& electrodeNames
        );

        //- Append one data row: time followed by all values
        static void writeRow
        (
            OFstream& os,
            scalar time,
            const List<scalar>& values
        );

        //- Load a torso surface STL from the 'torsoSurface' key in dict.
        //  Sets faceCentres (face barycentres) and vtkDir
        //  (runTime.path()/postProcessing/ECG, created automatically).
        //  Returns a null autoPtr if the key is absent.
        static autoPtr<triSurface> loadSurface
        (
            const dictionary& dict,
            const Time& runTime,
            pointField& faceCentres,
            fileName& vtkDir
        );

        //- Write face-centred ECG potentials as VTK ASCII POLYDATA.
        //  Filename: vtkDir/phiT_<timeName>.vtk
        //  fieldName is used for the SCALARS header (default: "phi").
        static void writeSurfaceVtk
        (
            const triSurface& surf,
            const List<scalar>& phi,
            const fileName& vtkDir,
            const word& timeName,
            const word& fieldName = "phi"
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
