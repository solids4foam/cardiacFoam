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

#include "conormalZeroFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{

conormalZeroFluxFvPatchScalarField::conormalZeroFluxFvPatchScalarField
(
    const fvPatch& p, const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    conductivityName_(),
    offsetFieldName_()
{
    gradient() = Zero;
}

conormalZeroFluxFvPatchScalarField::conormalZeroFluxFvPatchScalarField
(
    const fvPatch& p, const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    conductivityName_(dict.getOrDefault<word>("conductivity", word::null)),
    offsetFieldName_(dict.getOrDefault<word>("offsetField", word::null))
{
    gradient() = Zero;
    fvPatchScalarField::operator=(patchInternalField());
}

conormalZeroFluxFvPatchScalarField::conormalZeroFluxFvPatchScalarField
(
    const conormalZeroFluxFvPatchScalarField& ptf, const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    conductivityName_(ptf.conductivityName_),
    offsetFieldName_(ptf.offsetFieldName_)
{}

conormalZeroFluxFvPatchScalarField::conormalZeroFluxFvPatchScalarField
(
    const conormalZeroFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    conductivityName_(ptf.conductivityName_),
    offsetFieldName_(ptf.offsetFieldName_)
{}


void conormalZeroFluxFvPatchScalarField::setConductivity
(
    const word& conductivity,
    const word& offsetField
)
{
    if
    (
        (!conductivityName_.empty() && conductivityName_ != conductivity)
     || (!offsetFieldName_.empty() && offsetFieldName_ != offsetField)
    )
    {
        FatalErrorInFunction
            << "conormalZeroFlux patch " << patch().name() << " of field "
            << internalField().name() << " reads conductivity "
            << conductivityName_ << " offsetField " << offsetFieldName_
            << ", the solver uses conductivity " << conductivity
            << " offsetField " << offsetField
            << exit(FatalError);
    }
    conductivityName_ = conductivity;
    offsetFieldName_ = offsetField;
}


void conormalZeroFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (conductivityName_.empty())
    {
        FatalErrorInFunction
            << "conormalZeroFlux patch " << patch().name() << " of field "
            << internalField().name() << " has no conductivity"
            << exit(FatalError);
    }

    vectorField gradP
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + internalField().name() + ")"
        ).patchInternalField()
    );
    if (!offsetFieldName_.empty())
    {
        gradP +=
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + offsetFieldName_ + ")"
            ).patchInternalField();
    }
    const tensorField G
    (
        patch().lookupPatchField<volTensorField, tensor>(conductivityName_)
       .patchInternalField()
    );
    const vectorField n(patch().nf());
    const vectorField Gn(G & n);
    const scalarField nGn(n & Gn);
    if (min(nGn) <= 0)
    {
        FatalErrorInFunction
            << "conormalZeroFlux patch " << patch().name() << " of field "
            << internalField().name() << ": n.G.n <= 0 for conductivity "
            << conductivityName_
            << exit(FatalError);
    }

    gradient() = -((Gn - nGn*n) & gradP)/nGn;
    if (!offsetFieldName_.empty())
    {
        gradient() -=
            patch().lookupPatchField<volScalarField, scalar>(offsetFieldName_)
           .snGrad();
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void conormalZeroFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    if (!conductivityName_.empty())
    {
        os.writeEntry("conductivity", conductivityName_);
    }
    if (!offsetFieldName_.empty())
    {
        os.writeEntry("offsetField", offsetFieldName_);
    }
    fvPatchScalarField::writeValueEntry(os);
}


wordList conormalWallPatchTypes(const fvMesh& mesh, const word& trace)
{
    if (trace != "zeroGradient" && trace != "conormal")
    {
        FatalErrorInFunction
            << "sealedWallTrace must be zeroGradient or conormal"
            << exit(FatalError);
    }
    const word physicalType
    (
        trace == "conormal"
      ? conormalZeroFluxFvPatchScalarField::typeName
      : word("zeroGradient")
    );
    wordList types(mesh.boundary().size());
    forAll(mesh.boundary(), patchI)
    {
        const word& t = mesh.boundary()[patchI].type();
        types[patchI] = polyPatch::constraintType(t) ? t : physicalType;
    }
    return types;
}


void setConormalWallConductivity
(
    volScalarField& field,
    const word& conductivity,
    const word& offsetField
)
{
    volScalarField::Boundary& bf = field.boundaryFieldRef();
    forAll(bf, patchI)
    {
        if (isA<conormalZeroFluxFvPatchScalarField>(bf[patchI]))
        {
            refCast<conormalZeroFluxFvPatchScalarField>(bf[patchI])
               .setConductivity(conductivity, offsetField);
        }
    }
}


makePatchTypeField(fvPatchScalarField, conormalZeroFluxFvPatchScalarField);

} // End namespace Foam

// ************************************************************************* //
