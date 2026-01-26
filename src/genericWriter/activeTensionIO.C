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
    if (!userList.size())
    {
        return defaultNames;
    }

    wordList filtered;
    forAll(userList, i)
    {
        if (userList[i] == "Ta")
        {
            filtered.append(userList[i]);
        }
        else
        {
            WarningInFunction
                << "Ignoring unknown active tension field '" << userList[i]
                << "'. Only 'Ta' is supported." << nl;
        }
    }

    return filtered.size() ? filtered : defaultNames;
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

        outFields[i].primitiveFieldRef() = Ta;
        outFields[i].correctBoundaryConditions();
    }
}

} // End namespace Foam
