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
    conormalZeroFluxFvPatchScalarField

Description
    Sealed-wall patch with the conormal-consistent boundary value.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef conormalZeroFluxFvPatchScalarField_H
#define conormalZeroFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class fvMesh;

class conormalZeroFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    word conductivityName_;
    word offsetFieldName_;

public:

    TypeName("conormalZeroFlux");

    conormalZeroFluxFvPatchScalarField
    (
        const fvPatch&, const DimensionedField<scalar, volMesh>&
    );
    conormalZeroFluxFvPatchScalarField
    (
        const fvPatch&, const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );
    conormalZeroFluxFvPatchScalarField
    (
        const conormalZeroFluxFvPatchScalarField&, const fvPatch&,
        const DimensionedField<scalar, volMesh>&, const fvPatchFieldMapper&
    );
    conormalZeroFluxFvPatchScalarField
    (
        const conormalZeroFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new conormalZeroFluxFvPatchScalarField(*this)
        );
    }
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new conormalZeroFluxFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

wordList conormalWallPatchTypes(const fvMesh& mesh, const word& trace);

} // End namespace Foam

#endif

// ************************************************************************* //
