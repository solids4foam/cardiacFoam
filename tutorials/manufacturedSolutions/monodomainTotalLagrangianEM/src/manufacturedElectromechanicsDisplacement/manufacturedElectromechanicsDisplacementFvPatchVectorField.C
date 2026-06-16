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
#include "IOdictionary.H"

namespace Foam
{

// Read the manufactured displacement amplitude from the (top-level)
// electroMechanicalProperties dictionary. This MUST be done at construction:
// the same read works during case set-up but, in a decomposed parallel run,
// fails at run time because the file handler then looks for the dictionary in
// the region-local processor*/constant directory where it does not exist. The
// amplitude is constant in time, so it is read once and cached.
static vector readManufacturedAmplitude(const objectRegistry& db)
{
    IOdictionary emProps
    (
        IOobject
        (
            "electroMechanicalProperties",
            db.time().constant(),
            db.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word emModel(emProps.get<word>("electroMechanicalModel"));
    const dictionary& emCoeffs = emProps.subDict(emModel + "Coeffs");

    return emCoeffs.subDict("electromechanicalVerificationModel")
       .get<vector>("amplitude");
}

manufacturedElectromechanicsDisplacementFvPatchVectorField::
manufacturedElectromechanicsDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    amplitude_(vector::zero)
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
    amplitude_(readManufacturedAmplitude(p.boundaryMesh().mesh()))
{
    Info<< "Creating " << type() << " boundary condition, amplitude = "
        << amplitude_ << endl;
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

    // amplitude_ is read once at construction (see readManufacturedAmplitude);
    // it is constant in time and must not be re-read here, since a run-time
    // dictionary read is not parallel-safe.
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
    // amplitude is deduced from electroMechanicalProperties, not part of the
    // patch's serialized state, so it is not written here.
    fixedDisplacementFvPatchVectorField::write(os);
}


makePatchTypeField
(
    fvPatchVectorField,
    manufacturedElectromechanicsDisplacementFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
