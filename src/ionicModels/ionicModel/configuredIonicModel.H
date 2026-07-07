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

    virtual void configureApexBaseBandsHeterogeneity
    (
        const scalarField& apexDist,
        const dictionary& dict
    ) override
    {
        ionicHeterogeneityOrchestrator::configureApexBaseBandsHeterogeneityImpl
        (
            *this, apexDist, dict, HETEROGENEOUS_CONSTANTS_
        );
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
