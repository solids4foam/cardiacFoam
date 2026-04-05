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
    Foam::BuenoOrovioGPU

Description
    GPU-friendly batched version of the reduced Bueno-Orovio 2008 ionic model.

SourceFiles
    BuenoOrovioGPU.C
\*---------------------------------------------------------------------------*/

#ifndef BuenoOrovioGPU_H
#define BuenoOrovioGPU_H

#include "../ionicModel/configuredBatchedIonicModel.H"
#include "BuenoOrovio_2008Names.H"
#include "OFstream.H"
#include "volFields.H"
#include "overrideTypeName.H"

namespace Foam
{

class BuenoOrovioGPU
:
    public configuredBatchedIonicModel
{
    // Private Member Functions

        BuenoOrovioGPU(const BuenoOrovioGPU&) = delete;
        void operator=(const BuenoOrovioGPU&) = delete;


protected:
        virtual scalar ionicCurrentFromEvaluation
        (
            const scalarUList& algebraicValues
        ) const override
        {
            return algebraicValues[Jion] * 85.7;
        }

        virtual void evaluateState
        (
            const scalar modelTime,
            const scalarUList& stateValues,
            scalarUList& rateValues,
            scalarUList& algebraicValues
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


public:

        OverrideTypeName("BuenoOrovioGPU");

        virtual List<word> supportedTissueTypes() const override;
        virtual List<word> supportedDimensions() const override
        {
            return {};
        }

        BuenoOrovioGPU
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        );

        virtual ~BuenoOrovioGPU();

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
