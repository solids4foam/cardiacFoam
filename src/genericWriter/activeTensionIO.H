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

#ifndef activeTensionIO_H
#define activeTensionIO_H

#include "OFstream.H"
#include "PtrList.H"
#include "scalarField.H"
#include "wordList.H"
#include "volFields.H"
#include "error.H"

namespace Foam
{

class activeTensionIO
{
public:

    enum ColumnSource
    {
        COL_STATE     = 0,
        COL_ALGEBRAIC = 1,
        COL_RATE      = 2
    };

    struct SelectionPlan
    {
        wordList names;
        List<label> source;
        List<label> index;
    };

    struct SelectedMapCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr   = nullptr;
        label nStates = -1;
        label nAlg    = -1;
        wordList exportedNames;
        SelectionPlan plan;
    };

    struct FullPlanCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr   = nullptr;
        label nStates = -1;
        label nAlg    = -1;
        SelectionPlan plan;
    };

    struct ExportedNamesCache
    {
        bool valid = false;
        const void* stateNamesPtr = nullptr;
        const void* algNamesPtr   = nullptr;
        label nStates = -1;
        label nAlg    = -1;
        wordList requested;
        wordList filtered;
    };

    // Write column header for all states + algebraic + rates
    static void writeHeader
    (
        OFstream& os,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        FullPlanCache& fullPlanCache
    );

    // Write column header for a user-selected subset
    static void writeSelectedHeader
    (
        OFstream& os,
        const wordList& exportedNames
    );

    // Write one row: all states, algebraic, rates (integration point 0)
    static void write
    (
        scalar t,
        OFstream& os,
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES
    );

    // Write one row for a user-selected subset (integration point 0)
    static void writeSelected
    (
        const scalar t,
        OFstream& os,
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const wordList& exportedNames,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        SelectedMapCache& selectedPlanCache,
        const PtrList<scalarField>& RATES
    );

    // Return the filtered (valid) subset of userList against available names
    static const wordList& exportedFieldNamesRef
    (
        const wordList& userList,
        const char* const stateNames[],
        label nStates,
        const char* const algNames[],
        label nAlg,
        ExportedNamesCache& cache
    );

    // Export selected fields to volScalarField output list
    static void exportStateFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& exportedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlanCache,
        PtrList<volScalarField>& outFields
    );

    // Print selected fields for a given cell to Info
    static void debugPrintFields
    (
        const PtrList<scalarField>& STATES,
        const PtrList<scalarField>& ALGEBRAIC,
        const PtrList<scalarField>& RATES,
        const wordList& printedNames,
        const char* const stateNames[],
        int nStates,
        const char* const algNames[],
        int nAlg,
        SelectedMapCache& selectedPlanCache,
        label cellI,
        scalar t1,
        scalar t2,
        scalar step
    );
};

} // End namespace Foam

#endif
