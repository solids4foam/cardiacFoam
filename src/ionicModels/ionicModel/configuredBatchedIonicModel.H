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
    Foam::configuredBatchedIonicModel

Description
    Shared metadata/configuration layer for batched GPU-oriented ionic models.
    It keeps the per-model class focused on model math by centralizing
    constants storage, metadata hooks, and affine Vm-state mapping.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef configuredBatchedIonicModel_H
#define configuredBatchedIonicModel_H

#include "batchedIonicModel.H"
#include "ionicModelFamilyInfo.H"
#include "error.H"

namespace Foam
{

class configuredBatchedIonicModel
:
    public batchedIonicModel
{
    template<class Model>
    friend class BatchedKernelModelHooksView;

protected:

    mutable scalarField CONSTANTS_;
    mutable PtrList<scalarField> HETEROGENEOUS_CONSTANTS_;
    ionicModelFamilyInfo metadata_;

    configuredBatchedIonicModel
    (
        const dictionary& dict,
        const label nIntegrationPoints,
        const scalar initialDeltaT,
        const Switch solveVmWithinODESolver,
        const label nStates,
        const label nAlgebraics,
        const ionicModelFamilyInfo& metadata
    )
    :
        batchedIonicModel
        (
            dict,
            nIntegrationPoints,
            initialDeltaT,
            solveVmWithinODESolver,
            nStates,
            nAlgebraics
        ),
        CONSTANTS_(metadata.nConstants, 0.0),
        HETEROGENEOUS_CONSTANTS_(),
        metadata_(metadata)
    {}

    bool hasHeterogeneousConstants() const
    {
        return !HETEROGENEOUS_CONSTANTS_.empty();
    }

    scalarField& constants(const label integrationPtI) const
    {
        if (hasHeterogeneousConstants())
        {
            if
            (
                integrationPtI < 0
             || integrationPtI >= HETEROGENEOUS_CONSTANTS_.size()
            )
            {
                FatalErrorInFunction
                    << "Invalid heterogeneous constants integration point "
                    << integrationPtI << " for ionic model " << type()
                    << ". Valid range is [0, "
                    << HETEROGENEOUS_CONSTANTS_.size() - 1 << "]."
                    << exit(FatalError);
            }

            return HETEROGENEOUS_CONSTANTS_[integrationPtI];
        }

        return CONSTANTS_;
    }

    void flattenHeterogeneousConstants(scalarField& flattened) const
    {
        if (!hasHeterogeneousConstants())
        {
            flattened.clear();
            return;
        }

        flattened.setSize
        (
            HETEROGENEOUS_CONSTANTS_.size()*metadata_.nConstants
        );

        forAll(HETEROGENEOUS_CONSTANTS_, cellI)
        {
            const scalarField& cellConstants = HETEROGENEOUS_CONSTANTS_[cellI];

            if (cellConstants.size() != metadata_.nConstants)
            {
                FatalErrorInFunction
                    << "Invalid heterogeneous constants size "
                    << cellConstants.size() << " for ionic model " << type()
                    << " at integration point " << cellI
                    << ". Expected " << metadata_.nConstants << "."
                    << exit(FatalError);
            }

            for (label constantI = 0; constantI < metadata_.nConstants; ++constantI)
            {
                flattened[cellI*metadata_.nConstants + constantI] =
                    cellConstants[constantI];
            }
        }
    }

public:

    virtual void configureIonicHeterogeneity
    (
        const scalarField& transmuralDistance,
        const dictionary& heterogeneityDict
    ) override
    {
        const wordList tissues = supportedTissueTypes();
        const bool supportsEndoMCellEpi =
            tissues.found("endocardialCells")
         && tissues.found("mCells")
         && tissues.found("epicardialCells");

        if (!supportsEndoMCellEpi)
        {
            FatalErrorInFunction
                << "ionicHeterogeneity was requested for batched ionic model "
                << type()
                << ", but this model does not support endocardial/M/epicardial "
                << "transmural heterogeneity. Supported batched heterogeneous "
                << "models are BuenoOroviocompactBatched, TNNPcompactBatched, "
                << "ToRORd_dynClcompactBatched, and TWorldcompactBatched. "
                << type() << " supports tissue types " << tissues << "."
                << exit(FatalError);
        }

        configureTransmuralBandHeterogeneity
        (
            transmuralDistance,
            heterogeneityDict,
            HETEROGENEOUS_CONSTANTS_
        );
    }

protected:

    virtual ionicModelIO::VmTransform ioVmTransform() const override
    {
        return metadata_.vmTransform;
    }

    virtual const scalarField* ioConstantsPtr() const override
    {
        return &CONSTANTS_;
    }

    virtual scalarField* ioMutableConstantsPtr() const override
    {
        return &CONSTANTS_;
    }

    virtual label ioNumConstants() const override
    {
        return metadata_.nConstants;
    }

    virtual const char* const* ioConstantNames() const override
    {
        return metadata_.constantNames;
    }

    virtual const char* const* ioStateNames() const override
    {
        return metadata_.stateNames;
    }

    virtual const char* const* ioAlgebraicNames() const override
    {
        return metadata_.algebraicNames;
    }

    virtual label voltageStateIndex() const override
    {
        return metadata_.voltageStateI;
    }

    virtual scalar timeScaleFactor() const override
    {
        return metadata_.timeScale;
    }

    virtual scalar vmToState(const scalar vm) const override
    {
        return metadata_.vmInputScale*vm + metadata_.vmInputOffset;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
