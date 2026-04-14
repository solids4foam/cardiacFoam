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

#ifndef ElectromechanicalSignalProvider_H
#define ElectromechanicalSignalProvider_H

#include "scalar.H"
#include "label.H"

namespace Foam
{

//- Signals that an electromechanical provider can supply.
//  Can be expanded for mechanobiology feedback in the future.
enum class CouplingSignal
{
    VM,     //- Membrane potential (mV)
    CAI,    //- Intracellular Ca concentration (mM)
};

/*---------------------------------------------------------------------------*\
              Class ElectromechanicalSignalProvider Declaration
\*---------------------------------------------------------------------------*/

class ElectromechanicalSignalProvider
{
public:

    //- Destructor
    virtual ~ElectromechanicalSignalProvider() = default;

    // Member Functions

        //- Does this provider support the requested signal?
        virtual bool hasSignal(const CouplingSignal s) const = 0;

        //- Return the signal value at integration point i
        virtual scalar signal
        (
            const label i,
            const CouplingSignal s
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
