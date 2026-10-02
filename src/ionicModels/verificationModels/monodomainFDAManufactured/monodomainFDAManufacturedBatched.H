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

#ifndef monodomainFDAManufacturedBatched_H
#define monodomainFDAManufacturedBatched_H

#include "configuredBatchedIonicModel.H"
#include "monodomainFDAManufactured_2014Names.H"
#include "monodomainFDAManufacturedBatch.H"
#include "overrideTypeName.H"
#include "StimulusProtocolPOD.H"
#include "volFields.H"
#include <cstddef>

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class monodomainFDAManufacturedBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

    monodomainFDAManufacturedBatched
    (
        const monodomainFDAManufacturedBatched&
    ) = delete;
    void operator=(const monodomainFDAManufacturedBatched&) = delete;

    const scalarField* manufacturedSourceTermPtr_;
    scalar manufacturedSourceChi_;
    scalar manufacturedSourceCm_;
    label manufacturedSourceStart_;

#ifdef HAS_CUDA
    bool useDevice_;
    mutable batchedIonicCoreCUDA cuda_;
    double* dVmStateStart_;
    double* dVmStateRate_;
    std::size_t vmBufferSize_;
    scalarField vmStateStartHost_;
    scalarField vmStateRateHost_;
#endif

    scalar manufacturedSourceCorrection(const label cellI) const;

protected:

    virtual bool usesCompactSupport() const override final
    {
        return true;
    }

    virtual scalar ionicCurrentFromEvaluation
    (
        const scalarUList& algebraicValues
    ) const override final
    {
        return algebraicValues[Iion]/CONSTANTS_[Cm];
    }

    virtual scalar ionicCurrentFromHotPathSupport
    (
        const scalarUList& supportValues
    ) const override final
    {
        return supportValues[MONODOMAIN_FDA_BATCH_SUPPORT_Iion]/CONSTANTS_[Cm];
    }

    virtual List<word> supportedTissueTypes() const override final
    {
        return {};
    }

    virtual List<word> supportedDimensions() const override final
    {
        return {"1D", "2D", "3D"};
    }

    virtual scalarField constantsForTissue(const label dimension) const
        override final;

    virtual scalarField initialStatesForTissue(const label dimension) const
        override final;

    virtual void evaluateState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& algebraicValues
    ) const override final;

    virtual void evaluateHotPathState
    (
        const scalar modelTime,
        const scalarUList& stateValues,
        scalarUList& rateValues,
        scalarUList& supportValues
    ) const override final;

public:

    OverrideTypeName("monodomainFDAManufacturedBatched");

    monodomainFDAManufacturedBatched
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver = false
    );

    virtual ~monodomainFDAManufacturedBatched();

    virtual word tissueName() const override
    {
        return (tissue() == 1) ? "1D"
             : (tissue() == 2) ? "2D"
             : (tissue() == 3) ? "3D"
             :                   "unknown";
    }

    virtual word verificationFamily() const override
    {
        return "monodomainFDAManufactured";
    }

    virtual label geometricDimension() const override
    {
        return tissue();
    }

    virtual bool supportsVmRateCoupling() const override
    {
        return true;
    }

    virtual bool supportsIonicCurrentEvaluation() const override
    {
        return true;
    }

    virtual label nEqns() const override
    {
        return NUM_STATES;
    }

    virtual void solveODE
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    virtual void evaluateIonicCurrent
    (
        const scalar t,
        const scalarField& Vm,
        scalarField& Im
    ) override;

    virtual void setManufacturedSourceTerm
    (
        const scalarField& sourceTerm,
        scalar chi,
        scalar Cm,
        label sourceStart = 0
    ) override;

    virtual void clearManufacturedSourceTerm() override;

    virtual void prepareIOAccess
    (
        const wordList& requestedNames,
        const bool needsAlgebraics
    ) const override;

    virtual void importFields
    (
        const volScalarField& Vm,
        const wordList& fieldNames,
        const PtrList<volScalarField>& inFields
    ) override;

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        scalarField& dydt
    ) const override;

#ifdef HAS_CUDA
    void solveOnDevice
    (
        const scalar stepStartTime,
        const scalar deltaT,
        const scalarField& Vm,
        scalarField& Im
    );

    void setDeviceVmExtrapolant
    (
        const double modelTimeOffset,
        const label nCells
    );
#endif
};

} // End namespace Foam

#endif
