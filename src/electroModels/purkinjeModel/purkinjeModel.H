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
    Foam::purkinjeModel

Description
    Abstract base class for Purkinje network models.

    Follows the same runtime-selection pattern as ecgModel / ionicModel.
    Concrete subclasses implement evolve() to advance the 1-D Purkinje network
    one timestep and inject coupling currents at Purkinje-ventricular junction
    (PVJ) cell indices into externalStimulusCurrent.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    purkinjeModel.C
\*---------------------------------------------------------------------------*/

#ifndef purkinjeModel_H
#define purkinjeModel_H

#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class purkinjeModel Declaration
\*---------------------------------------------------------------------------*/

class purkinjeModel
{
protected:

    //- Transmembrane voltage owned by the parent electroModel
    const volScalarField& Vm_;

    //- Heart mesh
    const fvMesh& mesh_;


public:

    TypeName("purkinjeModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        purkinjeModel,
        dictionary,
        (
            const volScalarField& Vm,
            const dictionary&     dict,
            const scalar          initialDeltaT
        ),
        (Vm, dict, initialDeltaT)
    );

    // Constructors

        //- Construct from Vm field and coefficients dictionary
        purkinjeModel
        (
            const volScalarField& Vm,
            const dictionary&     dict
        );

    // Selectors

        //- Return a new purkinjeModel selected from dict entry 'purkinjeModel'
        static autoPtr<purkinjeModel> New
        (
            const volScalarField& Vm,
            const dictionary&     dict,
            const scalar          initialDeltaT
        );

    // Destructor
    virtual ~purkinjeModel() = default;

    // Member Functions

        //- Advance the 1D network one timestep and add coupling currents
        //  into externalStimulusCurrent at PVJ cell indices.
        virtual void evolve
        (
            scalar                t0,
            scalar                dt,
            const volScalarField& Vm,
            volScalarField&       externalStimulusCurrent
        ) = 0;

        //- Finalize any model-specific reporting at the end of the run
        virtual void end() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
