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


void conormalZeroFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const word& fieldName = internalField().name();
    const word gradName("grad(" + fieldName + ")");
    const word Gname
    (
        !conductivityName_.empty() ? conductivityName_
      : fieldName == "phiE"
     && db().foundObject<volTensorField>("conductivityExtracellular")
      ? word("conductivityExtracellular")
      : db().foundObject<volTensorField>("conductivityIntracellular")
      ? word("conductivityIntracellular")
      : word("conductivity")
    );
    const word offsetName
    (
        !offsetFieldName_.empty() ? offsetFieldName_
      : fieldName != "phiE"
     && Gname == "conductivityIntracellular"
     && db().foundObject<volScalarField>("phiE")
      ? word("phiE")
      : word::null
    );

    if
    (
        db().foundObject<volVectorField>(gradName)
     && db().foundObject<volTensorField>(Gname)
     && (
            offsetName.empty()
         || db().foundObject<volVectorField>("grad(" + offsetName + ")")
        )
    )
    {
        vectorField gradP
        (
            patch().lookupPatchField<volVectorField, vector>(gradName)
           .patchInternalField()
        );
        if (!offsetName.empty())
        {
            gradP +=
                patch().lookupPatchField<volVectorField, vector>
                (
                    "grad(" + offsetName + ")"
                ).patchInternalField();
        }
        const tensorField G
        (
            patch().lookupPatchField<volTensorField, tensor>(Gname)
           .patchInternalField()
        );
        const vectorField n(patch().nf());
        const vectorField Gn(G & n);
        const scalarField nGn(n & Gn);

        gradient() = -((Gn - nGn*n) & gradP)/max(nGn, SMALL);
        if (!offsetName.empty())
        {
            gradient() -=
                patch().lookupPatchField<volScalarField, scalar>(offsetName)
               .snGrad();
        }
    }
    else
    {
        gradient() = Zero;
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


makePatchTypeField(fvPatchScalarField, conormalZeroFluxFvPatchScalarField);

} // End namespace Foam

// ************************************************************************* //
