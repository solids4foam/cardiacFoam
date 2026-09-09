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
    Foam::configuredIonicModel

Description
    Shared heterogeneity configuration layer for non-batched scalar ionic
    models. Stores per-integration-point heterogeneous constants and
    forwards heterogeneity setup to ionicHeterogeneityOrchestrator.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef configuredIonicModel_H
#define configuredIonicModel_H

#include "ionicModel.H"
#include "ionicHeterogeneityOrchestrator.H"
#include "restartStateIO.H"

#include <fstream>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class configuredIonicModel Declaration
\*---------------------------------------------------------------------------*/

class configuredIonicModel
:
    public ionicModel
{
protected:

    //- Per-integration-point constants when transmural/named-region/
    //  cellZone or apex-base heterogeneity is configured
    mutable PtrList<scalarField> HETEROGENEOUS_CONSTANTS_;

public:

    //- Inherit ionicModel's constructor unchanged
    using ionicModel::ionicModel;

    virtual bool supportsRestartState() const override
    {
        return ioStatesPtr() != nullptr;
    }

    virtual bool readRestartState(const fvMesh& mesh) override
    {
        PtrList<scalarField>* states =
            const_cast<PtrList<scalarField>*>(ioStatesPtr());
        if (states == nullptr || states->empty())
        {
            return false;
        }

        const fileName statePath = restartStateIO::path(mesh, type() + "State");
        if (!isFile(statePath))
        {
            return false;
        }

        std::ifstream is(statePath.c_str(), std::ios::binary);
        if (!is)
        {
            FatalErrorInFunction
                << "Cannot read restart state file " << statePath
                << exit(FatalError);
        }
        restartStateIO::validateHeader
        (
            type(), ioNumStates(), states->size(), is, statePath
        );
        forAll(*states, integrationPtI)
        {
            forAll((*states)[integrationPtI], stateI)
            {
                (*states)[integrationPtI][stateI] =
                    restartStateIO::readScalar(is, statePath);
                restartStateIO::checkValue
                (
                    (*states)[integrationPtI][stateI], statePath
                );
            }
        }
        return true;
    }

    virtual void writeRestartState(const fvMesh& mesh) const override
    {
        const PtrList<scalarField>* states = ioStatesPtr();
        if (states == nullptr || states->empty())
        {
            return;
        }

        const fileName statePath = restartStateIO::path(mesh, type() + "State");
        std::ofstream os(statePath.c_str(), std::ios::binary | std::ios::trunc);
        if (!os)
        {
            FatalErrorInFunction
                << "Cannot write restart state file " << statePath
                << exit(FatalError);
        }
        restartStateIO::writeHeader
        (
            type(), ioNumStates(), states->size(), os
        );
        forAll(*states, integrationPtI)
        {
            forAll((*states)[integrationPtI], stateI)
            {
                restartStateIO::checkValue
                (
                    (*states)[integrationPtI][stateI], statePath
                );
                restartStateIO::writeScalar
                (
                    os, (*states)[integrationPtI][stateI], statePath
                );
            }
        }
    }

    virtual void refreshRestartState(const fvMesh& mesh) override
    {
        const volScalarField& Vm =
            mesh.lookupObject<volScalarField>("Vm");
        scalarField Im(Vm.size(), 0.0);
        evaluateIonicCurrent(mesh.time().value(), Vm, Im);
    }

    virtual void configureIonicHeterogeneity
    (
        const scalarField& transmuralDistance,
        const dictionary& heterogeneityDict
    ) override
    {
        ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity
        (
            *this, transmuralDistance, heterogeneityDict, HETEROGENEOUS_CONSTANTS_
        );
    }

    virtual void configureGradientAxisHeterogeneity
    (
        const scalarField& fieldValues,
        const dictionary& dict
    ) override
    {
        ionicHeterogeneityOrchestrator::configureGradientAxisHeterogeneity
        (
            *this, fieldValues, dict, HETEROGENEOUS_CONSTANTS_
        );
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
