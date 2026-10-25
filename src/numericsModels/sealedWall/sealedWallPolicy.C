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

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "sealedWallPolicy.H"
#include "labelIOList.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

bool Foam::sealedWallPolicy::isSealed(const fvPatchScalarField& psiP)
{
    return
        !psiP.coupled() && isA<zeroGradientFvPatchScalarField>(psiP);
}


Foam::labelList Foam::sealedWallPolicy::wallStates
(
    const volScalarField& psi,
    const bool seal
)
{
    labelList result(psi.boundaryField().size(), label(stock));
    if (!seal)
    {
        return result;
    }
    forAll(result, patchi)
    {
        const fvPatchScalarField& psiP = psi.boundaryField()[patchi];
        if (isSealed(psiP))
        {
            result[patchi] = sealed;
        }
        else if (!psiP.coupled() && isA<fixedValueFvPatchScalarField>(psiP))
        {
            result[patchi] = value;
        }
    }
    return result;
}


bool Foam::sealedWallPolicy::anyWall(const labelUList& states)
{
    for (const label s : states)
    {
        if (s != stock)
        {
            return true;
        }
    }
    return false;
}


Foam::word Foam::sealedWallPolicy::fieldKey(const word& gradName)
{
    std::string key(gradName);

    if (key.size() > 6 && key.compare(0, 5, "grad(") == 0 && key.back() == ')')
    {
        key = key.substr(5, key.size() - 6);
    }

    const auto dot = key.rfind(".component(");
    if (dot != std::string::npos && key.back() == ')')
    {
        const std::string index
        (
            key.substr(dot + 11, key.size() - dot - 12)
        );
        if
        (
            !index.empty()
         && index.find_first_not_of("0123456789") == std::string::npos
        )
        {
            key.erase(dot);
        }
    }

    return word(key, false);
}


Foam::word Foam::sealedWallPolicy::statesName(const word& key)
{
    return "wallGradientPatches_" + key;
}


void Foam::sealedWallPolicy::publish
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelUList& states
)
{
    if (states.size() != mesh.boundary().size())
    {
        FatalErrorInFunction
            << "Wall-gradient states of " << fieldName << " on mesh "
            << mesh.name() << " have " << states.size() << " entries for "
            << mesh.boundary().size() << " patches" << exit(FatalError);
    }

    const word name(statesName(fieldName));
    labelIOList* ptr = mesh.getObjectPtr<labelIOList>(name);
    if (ptr)
    {
        *ptr = states;
        return;
    }

    regIOobject::store
    (
        new labelIOList
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            states
        )
    );
}


void Foam::sealedWallPolicy::publish(const volScalarField& psi, const bool seal)
{
    publish(psi.mesh(), psi.name(), wallStates(psi, seal));
}


const Foam::labelList& Foam::sealedWallPolicy::states
(
    const fvMesh& mesh,
    const word& key
)
{
    const labelIOList* ptr = mesh.findObject<labelIOList>(statesName(key));
    if (!ptr)
    {
        FatalErrorInFunction
            << "No wall-gradient states for " << key << " on mesh "
            << mesh.name() << ". The owner of a field with sealed walls"
            << " publishes them; sealedWallLeastSquares applies only to"
            << " such fields." << exit(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    if (ptr->size() != patches.size())
    {
        FatalErrorInFunction
            << "Wall-gradient states of " << key << " on mesh " << mesh.name()
            << " have " << ptr->size() << " entries for " << patches.size()
            << " patches" << exit(FatalError);
    }

    forAll(*ptr, patchi)
    {
        const label s = (*ptr)[patchi];
        if (s < stock || s > value)
        {
            FatalErrorInFunction
                << "Wall-gradient state " << s << " of " << key
                << " on patch " << patches[patchi].name()
                << " is not 0 (stock), 1 (sealed) or 2 (value)"
                << exit(FatalError);
        }
        if
        (
            s != stock
         && (
                patches[patchi].coupled()
             || polyPatch::constraintType(patches[patchi].type())
            )
        )
        {
            FatalErrorInFunction
                << "Patch " << patches[patchi].name() << " ("
                << patches[patchi].type() << ") of mesh " << mesh.name()
                << " cannot be a sealed or value wall of " << key
                << exit(FatalError);
        }
    }

    return *ptr;
}


// ************************************************************************* //
