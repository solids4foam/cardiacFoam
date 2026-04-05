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
    Foam::ORdGPU

Description
    Batched GPU-oriented path for the ORd 2011 ionic model.

SourceFiles
    ORdGPU.C
\*---------------------------------------------------------------------------*/

#ifndef ORdGPU_H
#define ORdGPU_H

#include "../ionicModel/configuredBatchedIonicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "ORd_2011Names.H"
#include "../ORd/ORdCompactSupport.H"
#include "overrideTypeName.H"

namespace Foam
{

class ORdGPU
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

protected:

    enum class FullEvaluationBackend : unsigned char
    {
        generated,
        optimized
    };

    FullEvaluationBackend fullEvaluationBackend_;

    ORdGPU(const ORdGPU&) = delete;
    void operator=(const ORdGPU&) = delete;

    ORdGPU
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver,
        const FullEvaluationBackend fullEvaluationBackend
    );

    virtual label activeHotPathSupportSize
    (
        const bool useExplicitEuler
    ) const override
    {
        return useExplicitEuler ? 1 : nHotPathSupport();
    }

    virtual bool usesCompactSupport() const override
    {
        return true;
    }

    virtual scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const override
    {
        return algebraicValues[Iion_cm];
    }

    virtual scalar ionicCurrentFromHotPathSupport
    (
        const scalarUList& supportValues
    ) const override;

    virtual void stabilizeStateValues(scalarUList& stateValues) const override;

    virtual void stabilizeCellState(const label cellI) const override;

    virtual void evaluateState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const override;

    virtual void evaluateHotPathState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& supportValues
    ) const override;

    virtual bool rushLarsenParameters
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& algebraicValues,
        scalar& steadyState,
        scalar& tau
    ) const override;

    virtual bool rushLarsenParametersFromHotPathSupport
    (
        const label stateI,
        const scalarUList& stateValues,
        const scalarUList& rateValues,
        const scalarUList& supportValues,
        scalar& steadyState,
        scalar& tau
    ) const override;

public:

    OverrideTypeName("ORdGPU");

    virtual List<word> supportedTissueTypes() const override;
    virtual List<word> supportedDimensions() const override
    {
        return {};
    }

    ORdGPU
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~ORdGPU();

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override;

    virtual void sweepCurrent
    (
        const word& currentName,
        scalar Vmin,
        scalar Vmax,
        label nPts,
        const fileName& outputFile
    ) const override;

    virtual wordList availableSweepCurrents() const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
