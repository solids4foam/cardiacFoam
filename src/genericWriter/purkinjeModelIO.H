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
    purkinjeModelIO

Description
    Static utility class for Purkinje network model I/O operations.

Author
    Simao Nieto de Castro. All rights reserved.
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
#include "PtrList.H"
#include "globalIOFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class purkinjeModelIO Declaration
\*---------------------------------------------------------------------------*/

class purkinjeModelIO
{
public:

    struct ResolvedTokens
    {
        wordList networkTokens;
        wordList ionicTokens;
        wordList unknownTokens;
    };

    // Constructors

        //- Not instantiable – all methods are static
        purkinjeModelIO() = delete;


    // Static Member Functions

        static ResolvedTokens filterTokens
        (
            const wordList& userExport,
            const char* const stateNames[],
            int nStates,
            const char* const algNames[],
            int nAlg
        );

        //- Open/create a time-series output file and write the column header.
        //  Creates outDir if it does not exist. A run starting after time 0
        //  keeps the rows up to startTime that an earlier run wrote under the
        //  same header, so a restart continues the series.
        static autoPtr<OFstream> openTimeSeries
        (
            const fileName& outDir,
            const word& filename,
            const wordList& columnNames,
            scalar startTime
        );

        //- Entries up to startTime of an existing .vtk.series file
        static void readVTKSeries
        (
            const fileName& seriesPath,
            scalar startTime,
            scalarList& times,
            wordList& filenames
        );

        //- Write a graph field held whole on every rank into the current
        //  time directory of the undecomposed case, from the master.
        static void writeGlobalField(scalarGlobalIOField& field);

        //- Read a graph field written by writeGlobalField from the field's
        //  time directory (its instance, the start time for a field
        //  constructed at start-up) of the undecomposed case: the master
        //  reads it and every rank receives it. Returns false, leaving the
        //  field unchanged, when there is no such file.
        static bool readGlobalField(scalarGlobalIOField& field);

        //- Append one data row: time followed by all values
        static void writeRow
        (
            OFstream& os,
            scalar time,
            const List<scalar>& values
        );

        //- (Re)write a .vtk.series JSON sidecar so that ParaView (≥5.6)
        //  displays correct physical time values for legacy VTK files.
        //  Pass all accumulated (time, filename) pairs each call.
        static void writeVTKSeries
        (
            const fileName&   seriesPath,
            const scalarList& times,
            const wordList&   filenames
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
        //  includeVmIion : false for ionic-model-less (eikonal) 1D solvers.
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
            const scalarField&  terminalSource,
            bool                includeVmIion = true,
            const PtrList<scalarField>& ionicFields = PtrList<scalarField>(),
            const wordList&             ionicFieldNames = wordList()
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
