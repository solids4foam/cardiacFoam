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
#include "pointField.H"
#include "labelList.H"
#include "scalarField.H"
#include "scalarList.H"

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

        //- Write one VTK legacy PolyData file for the Purkinje graph at a
        //  given simulation time.  Each call produces an independent
        //  purkinjeNetwork_<timeStepIndex>.vtk file in outDir that ParaView
        //  automatically groups into a time series using zero-padded indices.
        //
        //  nodeLocations : 3-D position of every graph node (size = nNodes).
        //                  If empty the function is a no-op.
        //  edgeNodeA/B   : parallel lists giving the two endpoint node indices
        //                  of every edge (size = nEdges each).
        //  Vm1D          : membrane potential at every node (BuenoOrovio scale).
        //  Iion1D        : ionic current at every node.
        //  pvjNodes      : node indices of the PVJ terminals.
        //  terminalSource: volumetric source at each PVJ (size = pvjNodes.size()).
        //- (Re)write a .vtk.series JSON sidecar so that ParaView (≥5.6)
        //  displays correct physical time values for legacy VTK files.
        //  Pass all accumulated (time, filename) pairs each call.
        static void writeVTKSeries
        (
            const fileName&   seriesPath,
            const scalarList& times,
            const wordList&   filenames
        );

        static void writeVTK
        (
            const fileName&     outDir,
            const word&         timeName,
            const label         timeStep,
            const pointField&   nodeLocations,
            const labelList&    edgeNodeA,
            const labelList&    edgeNodeB,
            const scalarField&  Vm1D,
            const scalarField&  Iion1D,
            const labelList&    pvjNodes,
            const scalarField&  terminalSource
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
