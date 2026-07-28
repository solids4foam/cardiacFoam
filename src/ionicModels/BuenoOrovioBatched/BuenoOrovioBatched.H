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
    Foam::BuenoOrovioBatched

Description
    GPU-friendly batched version of the reduced Bueno-Orovio 2008 ionic model.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    BuenoOrovioBatched.C
\*---------------------------------------------------------------------------*/

#ifndef BuenoOrovioGPU_H
#define BuenoOrovioGPU_H

#include "configuredBatchedIonicModel.H"
#include "BuenoOrovio_2008Names.H"
#include "BuenoOrovio_2008Batch.H"
#include "StimulusProtocolPOD.H"
#include "OFstream.H"
#include "volFields.H"
#include "overrideTypeName.H"

#ifdef HAS_CUDA
#include "batchedIonicCoreCUDA.H"
#endif

namespace Foam
{

class BuenoOrovioBatched
:
    public configuredBatchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;


        BuenoOrovioBatched(const BuenoOrovioBatched&) = delete;
        void operator=(const BuenoOrovioBatched&) = delete;

protected:

        mutable StimulusProtocolPOD stimulusPOD_;

#ifdef HAS_CUDA
        bool useDevice_;
        mutable batchedIonicCoreCUDA cuda_;
#endif // HAS_CUDA

        virtual scalar ionicCurrentFromEvaluation
        (
            const scalarUList& algebraicValues
        ) const override final
        {
            return algebraicValues[Jion] * 85.7;
        }

        virtual scalar ionicCurrentFromHotPathSupport
        (
            const scalarUList& supportValues
        ) const override final;

        virtual scalarField constantsForTissue(const label tissueFlag) const
            override final;

        virtual scalarField initialStatesForTissue(const label tissueFlag) const
            override final;

        virtual bool usesCompactSupport() const override final
        {
            return true;
        }

        virtual void evaluateState
        (
            const scalar modelTime,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& algebraicValues
        ) const override final;

        virtual void evaluateState
        (
            const label cellI,
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

        virtual void evaluateHotPathState
        (
            const label cellI,
            const scalar modelTime,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& supportValues
        ) const override final;

        virtual bool rushLarsenParametersFromHotPathSupport
        (
            const label stateI,
            const scalarUList& stateValues,
            const scalarUList& rateValues,
            const scalarUList& supportValues,
            scalar& steadyState,
            scalar& tau
        ) const override final;

        virtual bool rushLarsenParametersFromHotPathSupport
        (
            const label cellI,
            const label stateI,
            const scalarUList& stateValues,
            const scalarUList& rateValues,
            const scalarUList& supportValues,
            scalar& steadyState,
            scalar& tau
        ) const override final;


public:

        OverrideTypeName("BuenoOrovioBatched");

        virtual List<word> supportedTissueTypes() const override final;
        virtual List<word> supportedDimensions() const override final
        {
            return {};
        }

        BuenoOrovioBatched
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );

        virtual ~BuenoOrovioBatched();

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        ) override;

        virtual void prepareIOAccess
        (
            const wordList& requestedNames,
            const bool needsAlgebraics
        ) const override final;

        virtual void importFields
        (
            const volScalarField& Vm,
            const wordList& fieldNames,
            const PtrList<volScalarField>& inFields
        ) override final;

        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const override final;
    };

    class BuenoOroviocompactBatched final
    :
        public BuenoOrovioBatched
    {
        BuenoOroviocompactBatched(const BuenoOroviocompactBatched&) = delete;
        void operator=(const BuenoOroviocompactBatched&) = delete;

    public:

        OverrideTypeName("BuenoOroviocompactBatched");

        BuenoOroviocompactBatched
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );

        virtual ~BuenoOroviocompactBatched() = default;

        virtual void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        ) override;

    protected:

#ifdef HAS_CUDA
        void solveOnDevice
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im
        );
#endif
    };

} // End namespace Foam

#endif

// ************************************************************************* //
