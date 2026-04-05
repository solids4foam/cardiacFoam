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
\*---------------------------------------------------------------------------*/

#ifndef configuredBatchedIonicModel_H
#define configuredBatchedIonicModel_H

#include "batchedIonicModel.H"
#include "ionicModelFamilyInfo.H"

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
        metadata_(metadata)
    {}

    virtual ionicModelIO::VmTransform ioVmTransform() const override
    {
        return metadata_.vmTransform;
    }

    virtual const scalarField* ioConstantsPtr() const override
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
