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
    Foam::ManufacturedElectromechanics

Description
    Algebraic manufactured active tension model driven by Vm and fibre stretch.

SourceFiles
    ManufacturedElectromechanics.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef ManufacturedElectromechanics_H
#define ManufacturedElectromechanics_H

#include "activeTensionModel.H"
#include "overrideTypeName.H"

#include "PtrList.H"
#include "scalarField.H"

#include "ManufacturedElectromechanics_2026Names.H"

namespace Foam
{

class ManufacturedElectromechanics
:
    public activeTensionModel
{
    // Private data

        //- Per-integration-point active tension state storage
        PtrList<scalarField> STATES_;

        //- Per-integration-point algebraic variable storage
        PtrList<scalarField> ALGEBRAIC_;

        //- Per-integration-point rate storage
        PtrList<scalarField> RATES_;

        //- Shared constant storage
        mutable scalarField CONSTANTS_;


protected:

        virtual const PtrList<scalarField>* ioStatesPtr() const override
        {
            return &STATES_;
        }

        virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
        {
            return &ALGEBRAIC_;
        }

        virtual const PtrList<scalarField>* ioRatesPtr() const override
        {
            return &RATES_;
        }

        virtual const scalarField* ioConstantsPtr() const override
        {
            return &CONSTANTS_;
        }

        virtual label ioNumConstants() const override
        {
            return NUM_CONSTANTS;
        }

        virtual const char* const* ioConstantNames() const override;

        virtual const char* const* ioStateNames() const override;

        virtual label ioNumStates() const override
        {
            return NUM_STATES;
        }

        virtual const char* const* ioAlgebraicNames() const override;

        virtual label ioNumAlgebraic() const override
        {
            return NUM_ALGEBRAIC;
        }

        virtual CouplingSignal driveSignal() const override
        {
            return CouplingSignal::VM;
        }

        //- Manufactured solution is posed directly on the OpenFOAM clock.
        virtual scalar timeScaleFactor() const override
        {
            return 1.0;
        }

        virtual void solveAtPoint
        (
            label i,
            scalar driveVal,
            scalar lambda,
            scalar& Ta
        ) override;


public:

    OverrideTypeName("ManufacturedElectromechanics");

    // Constructors

        ManufacturedElectromechanics
        (
            const dictionary& dict,
            const label nIntegrationPoints
        );

    // Destructor

        virtual ~ManufacturedElectromechanics() = default;


    // Member Functions

        virtual Requirements requirements() const override;

        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
