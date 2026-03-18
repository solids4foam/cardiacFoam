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

#ifndef CouplingSignalProvider_H
#define CouplingSignalProvider_H

#include "scalar.H"
#include "label.H"

namespace Foam
{

enum class CouplingSignal // this can be expanded for the necessary eletromechanic feedback. For mechanobiology, see later. 
{
    Vm,        // membrane potential mV
    Cai,       // intracellular Ca  mM
};

class CouplingSignalProvider
{
public:
    virtual ~CouplingSignalProvider() = default;

    // Does this provider support the requested signal?
    virtual bool hasSignal(const CouplingSignal s) const = 0;
    virtual scalar signal(const label i, const CouplingSignal s) const = 0;
};

} // End namespace Foam

#endif
