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
    Foam::ORdGPUOpt

Description
    Legacy runtime-selectable ORd alias that keeps the historical
    `ORdGPUOpt` name while selecting the optimized full-evaluation backend
    inside the shared `ORdGPU` architecture.

SourceFiles
    ORdGPUOpt.C
\*---------------------------------------------------------------------------*/

#ifndef ORdGPUOpt_H
#define ORdGPUOpt_H

#include "ORdGPU.H"
#include "overrideTypeName.H"

namespace Foam
{

class ORdGPUOpt
:
    public ORdGPU
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    ORdGPUOpt(const ORdGPUOpt&) = delete;
    void operator=(const ORdGPUOpt&) = delete;

public:

    OverrideTypeName("ORdGPUOpt");

    ORdGPUOpt
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~ORdGPUOpt();

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
