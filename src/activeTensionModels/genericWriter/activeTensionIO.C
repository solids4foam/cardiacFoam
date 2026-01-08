/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cardiacFoam
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Simple I/O helpers for active tension models
\*---------------------------------------------------------------------------*/

#include "activeTensionIO.H"

namespace Foam
{

wordList activeTensionIO::exportedFieldNames
(
    const wordList& userList,
    const wordList& defaultNames
)
{
    return userList.size() ? userList : defaultNames;
}


void activeTensionIO::writeHeader
(
    OFstream& os,
    const wordList& names
)
{
    os << "# t";
    forAll(names, i)
    {
        os << " " << names[i];
    }
    os << nl;
}


void activeTensionIO::write
(
    const scalar t,
    OFstream& os,
    const scalarField& values
)
{
    os << t;
    forAll(values, i)
    {
        os << " " << values[i];
    }
    os << nl;
}


void activeTensionIO::exportStateFields
(
    const scalarField& Ta,
    const wordList& exportedNames,
    PtrList<volScalarField>& outFields
)
{
    if (exportedNames.size() != outFields.size())
    {
        FatalErrorInFunction
            << "exportedNames.size() != outFields.size()"
            << abort(FatalError);
    }

    forAll(exportedNames, i)
    {
        if (exportedNames[i] != "Ta")
        {
            FatalErrorInFunction
                << "Unknown active tension field name '" << exportedNames[i]
                << "'. Only 'Ta' is supported."
                << abort(FatalError);
        }

        if (outFields[i].internalField().size() != Ta.size())
        {
            FatalErrorInFunction
                << "Ta size does not match outFields internal size."
                << abort(FatalError);
        }

        outFields[i].internalField() = Ta;
        outFields[i].correctBoundaryConditions();
    }
}

} // End namespace Foam
