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

Function
    readEcgElectrodes

Description
    Electrode names and positions of an ECG domain entry: its
    electrodePositions, else the electrodePositions shared by all entries of
    ecgDomains. No electrodes is fatal. The single reader of the electrodes,
    used by ecgDomain and by the lead-field build (leadFieldConductor).

SourceFiles
    ecgElectrodes.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef ecgElectrodes_H
#define ecgElectrodes_H

#include "dictionary.H"
#include "wordList.H"
#include "vectorList.H"

namespace Foam
{

//- Electrode names and positions, in the order of the dictionary
struct ecgElectrodes
{
    wordList names;
    List<vector> positions;
};

//- The electrodes of the ECG domain entry ecgDict; sharedPositionsPtr is
//  the electrodePositions of ecgDomains, or nullptr
ecgElectrodes readEcgElectrodes
(
    const dictionary& ecgDict,
    const dictionary* sharedPositionsPtr
);

} // End namespace Foam

#endif

// ************************************************************************* //
