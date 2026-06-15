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
\*---------------------------------------------------------------------------*/

#include "manufacturedElectromechanicsDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "manufacturedElectromechanicsReference.H"

namespace Foam
{

manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    amplitude_(vector(0.02, 0.02, 0.02))
{}


manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const manufacturedElectromechanicsDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_)
{}


manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    amplitude_(dict.lookupOrDefault<vector>("amplitude", vector(0.02, 0.02, 0.02)))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}

manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const manufacturedElectromechanicsDisplacementFvPatchVectorField& ptf
)
:
    fixedDisplacementFvPatchVectorField(ptf),
    amplitude_(ptf.amplitude_)
{}


manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const manufacturedElectromechanicsDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(ptf, iF),
    amplitude_(ptf.amplitude_)
{}


void manufacturedElectromechanicsDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    computeManufacturedElectromechanicsD
    (
        totalDisp(),
        patch().Cf(),
        this->db().time().value(),
        amplitude_
    );

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


void manufacturedElectromechanicsDisplacementFvPatchVectorField::write
(
    Ostream& os
) const
{
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;

    fixedDisplacementFvPatchVectorField::write(os);
}


makePatchTypeField
(
    fvPatchVectorField,
    manufacturedElectromechanicsDisplacementFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
