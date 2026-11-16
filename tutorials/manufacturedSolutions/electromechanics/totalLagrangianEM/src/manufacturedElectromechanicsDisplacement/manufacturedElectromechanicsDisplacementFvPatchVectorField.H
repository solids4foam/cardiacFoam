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
    Foam::manufacturedElectromechanicsDisplacementFvPatchVectorField

Description
    Boundary condition that imposes the manufactured electromechanics
    displacement field on a solid patch while preserving the solids4foam
    fixed-displacement handling.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedElectromechanicsDisplacementFvPatchVectorField_H
#define manufacturedElectromechanicsDisplacementFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"

namespace Foam
{

class manufacturedElectromechanicsDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private data

        //- Manufactured displacement amplitude. Deduced from
        //  electroMechanicalProperties in updateCoeffs() so it stays in sync
        //  with the manufactured body force; not read from the patch dict.
        vector amplitude_;

public:

    TypeName("manufacturedElectromechanicsDisplacement");

    manufacturedElectromechanicsDisplacementFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    manufacturedElectromechanicsDisplacementFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    manufacturedElectromechanicsDisplacementFvPatchVectorField
    (
        const manufacturedElectromechanicsDisplacementFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    manufacturedElectromechanicsDisplacementFvPatchVectorField
    (
        const manufacturedElectromechanicsDisplacementFvPatchVectorField&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new manufacturedElectromechanicsDisplacementFvPatchVectorField(*this)
        );
    }

    manufacturedElectromechanicsDisplacementFvPatchVectorField
    (
        const manufacturedElectromechanicsDisplacementFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new manufacturedElectromechanicsDisplacementFvPatchVectorField
            (
                *this,
                iF
            )
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
