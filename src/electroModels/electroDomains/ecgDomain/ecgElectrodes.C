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

#include "ecgElectrodes.H"

namespace Foam
{

ecgElectrodes readEcgElectrodes
(
    const dictionary& ecgDict,
    const dictionary* sharedPositionsPtr
)
{
    const dictionary* eDictPtr = ecgDict.findDict("electrodePositions");

    if (!eDictPtr)
    {
        eDictPtr = sharedPositionsPtr;
    }

    if (!eDictPtr || eDictPtr->empty())
    {
        FatalIOErrorInFunction(ecgDict)
            << "ECG domain " << ecgDict.dictName() << " has no electrodes: "
            << "give electrodePositions in the domain or in ecgDomains."
            << exit(FatalIOError);
    }

    ecgElectrodes electrodes;
    electrodes.names = eDictPtr->toc();
    electrodes.positions.setSize(electrodes.names.size());

    forAll(electrodes.names, i)
    {
        electrodes.positions[i] = eDictPtr->get<vector>(electrodes.names[i]);
    }

    return electrodes;
}

} // End namespace Foam

// ************************************************************************* //
