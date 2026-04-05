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
    Foam::ORdGPUCompact

Description
    Legacy runtime-selectable ORd alias that keeps the historical
    `ORdGPUCompact` name for backwards compatibility. The hot-path/backend
    architecture now lives in `ORdGPU`; scheme-level behavior is selected in
    the case configuration rather than by maintaining a separate compact class
    implementation.
\*---------------------------------------------------------------------------*/

#ifndef ORdGPUCompact_H
#define ORdGPUCompact_H

#include "ORdGPUOpt.H"
#include "overrideTypeName.H"

namespace Foam
{

class ORdGPUCompact
:
    public ORdGPUOpt
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    ORdGPUCompact(const ORdGPUCompact&) = delete;
    void operator=(const ORdGPUCompact&) = delete;

public:

    OverrideTypeName("ORdGPUCompact");

    ORdGPUCompact
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~ORdGPUCompact();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
