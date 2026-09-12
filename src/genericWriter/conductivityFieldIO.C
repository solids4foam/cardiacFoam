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

#include "conductivityFieldIO.H"

#include "PstreamReduceOps.H"
#include "Switch.H"

namespace Foam
{

namespace
{

tmp<volTensorField> mapToSolverMesh
(
    const volTensorField& supportField,
    const fvMesh& solverMesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const word& internalName
)
{
    if (&solverMesh == &supportMesh)
    {
        return tmp<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    internalName,
                    solverMesh.time().timeName(),
                    solverMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                supportField
            )
        );
    }

    if (!meshSubsetPtr || !meshSubsetPtr->hasSubMesh())
    {
        FatalErrorInFunction
            << "Conductivity field '" << supportField.name()
            << "' was read on support mesh '" << supportMesh.name()
            << "', but solver mesh '" << solverMesh.name()
            << "' is different and no fvMeshSubset mapping was supplied."
            << exit(FatalError);
    }

    tmp<volTensorField> tmapped(meshSubsetPtr->interpolate(supportField));
    tmapped.ref().rename(internalName);
    return tmapped;
}


tmp<volTensorField> widenSymmetricField
(
    const volSymmTensorField& field,
    const word& internalName
)
{
    return tmp<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                internalName,
                field.time().timeName(),
                field.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            field & tensor(I)
        )
    );
}


// Find the time instance holding fieldName under supportMesh, or return an
// empty word if it genuinely cannot be found anywhere (not even
// "constant"). OpenFOAM.com added findInstance's constant_fallback argument
// in 2406; before that it always falls back to "constant" internally, so an
// empty result never gets returned and the fallback has to be detected by
// checking whether the field actually exists at the resolved instance.
word findFieldInstance(const fvMesh& supportMesh, const word& fieldName)
{
    word instance = supportMesh.time().findInstance
    (
        supportMesh.dbDir(),
        fieldName,
        IOobject::READ_IF_PRESENT,
        word::null
#if OPENFOAM >= 2406
        , false
#endif
    );

#if OPENFOAM < 2406
    if (!instance.empty())
    {
        IOobject io
        (
            fieldName,
            instance,
            supportMesh.dbDir(),
            supportMesh,
            IOobject::READ_IF_PRESENT
        );

        if (!io.typeHeaderOk<volTensorField>(false, false, false))
        {
            instance = word::null;
        }
    }
#endif

    return instance;
}


void validateLegacySymmetry(const volTensorField& field)
{
    scalar maxAsymmetry = 0.0;
    scalar maxMagnitude = 0.0;

    forAll(field, celli)
    {
        maxAsymmetry = max
        (
            maxAsymmetry,
            mag(field[celli] - field[celli].T())
        );
        maxMagnitude = max(maxMagnitude, mag(field[celli]));
    }

    reduce(maxAsymmetry, maxOp<scalar>());
    reduce(maxMagnitude, maxOp<scalar>());

    const scalar tolerance = 1e-12*max(scalar(1), maxMagnitude);
    if (maxAsymmetry > tolerance)
    {
        FatalErrorInFunction
            << "Legacy conductivity field '" << field.name()
            << "' is a volTensorField but is not symmetric. max|D-D.T| = "
            << maxAsymmetry << ", tolerance = " << tolerance << "."
            << exit(FatalError);
    }
}

} // End anonymous namespace


tmp<volTensorField> readConductivityField
(
    const fvMesh& solverMesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const dictionary& coefficients,
    const conductivityFieldSpec& spec
)
{
    word source("uniform");

    if (coefficients.found("conductivitySource"))
    {
        source = word(coefficients.lookup("conductivitySource"));
    }
    else
    {
        WarningInFunction
            << "No conductivitySource specified in " << coefficients.name()
            << "; using the legacy default 'uniform'. Add either "
            << "'conductivitySource uniform;' or 'conductivitySource field;'."
            << endl;
    }

    if (source == "uniform")
    {
        if (!coefficients.found(spec.dictionaryEntry))
        {
            FatalIOErrorInFunction(coefficients)
                << "conductivitySource is 'uniform' but entry '"
                << spec.dictionaryEntry << "' is missing."
                << exit(FatalIOError);
        }

        return tmp<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    spec.internalName,
                    solverMesh.time().timeName(),
                    solverMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solverMesh,
                dimensionedTensor
                (
                    dimensionedSymmTensor
                    (
                        spec.dictionaryEntry,
                        coefficients
                    ) & tensor(I)
                )
            )
        );
    }

    if (source != "field")
    {
        FatalIOErrorInFunction(coefficients)
            << "Unknown conductivitySource '" << source
            << "'. Valid values are 'field' and 'uniform'."
            << exit(FatalIOError);
    }

    word fieldName(spec.fieldName);
    word instance = findFieldInstance(supportMesh, fieldName);

    if (instance.empty() && spec.internalName != spec.fieldName)
    {
        const word legacyInstance =
            findFieldInstance(supportMesh, spec.internalName);

        if (!legacyInstance.empty())
        {
            fieldName = spec.internalName;
            instance = legacyInstance;
            WarningInFunction
                << "Reading deprecated conductivity field name '" << fieldName
                << "'. Regenerate it as canonical field '" << spec.fieldName
                << "'." << endl;
        }
    }

    if (instance.empty())
    {
        FatalIOErrorInFunction(coefficients)
            << "conductivitySource is 'field' but canonical field '"
            << spec.fieldName << "' was not found at or before time "
            << supportMesh.time().timeName() << "."
            << exit(FatalIOError);
    }

    IOobject fieldIO
    (
        fieldName,
        instance,
        supportMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (fieldIO.typeHeaderOk<volSymmTensorField>(true))
    {
        volSymmTensorField symmetricField
        (
            IOobject
            (
                fieldName,
                instance,
                supportMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            supportMesh
        );

        tmp<volTensorField> tsupport
        (
            widenSymmetricField(symmetricField, spec.internalName)
        );

        Info<< "Conductivity source: field " << instance << '/' << fieldName
            << " (volSymmTensorField) -> " << spec.internalName << nl << endl;

        return mapToSolverMesh
        (
            tsupport(),
            solverMesh,
            supportMesh,
            meshSubsetPtr,
            spec.internalName
        );
    }

    if (fieldIO.typeHeaderOk<volTensorField>(true))
    {
        volTensorField legacyField
        (
            IOobject
            (
                fieldName,
                instance,
                supportMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            supportMesh
        );

        validateLegacySymmetry(legacyField);

        WarningInFunction
            << "Reading legacy volTensorField conductivity '" << fieldName
            << "'. Canonical conductivity fields are volSymmTensorField."
            << endl;

        return mapToSolverMesh
        (
            legacyField,
            solverMesh,
            supportMesh,
            meshSubsetPtr,
            spec.internalName
        );
    }

    FatalIOErrorInFunction(coefficients)
        << "Conductivity field '" << instance << '/' << fieldName
        << "' has class '" << fieldIO.headerClassName()
        << "'; expected volSymmTensorField."
        << exit(FatalIOError);

    return tmp<volTensorField>(nullptr);
}

} // End namespace Foam

// ************************************************************************* //
