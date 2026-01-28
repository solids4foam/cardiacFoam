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

void activeTensionIO::debugPrintFields
(
    const wordList& printedNames,
    const wordList& availableNames,
    const scalarField& values,
    label cellI,
    scalar t1,
    scalar t2,
    scalar step
)
{
    if (printedNames.empty())
    {
        return;
    }

    Info<< "DEBUG cell=" << cellI
        << " t=" << t1 << "->" << t2;

    if (step >= 0)
    {
        Info<< " step=" << step;
    }

    forAll(printedNames, k)
    {
        const word& name = printedNames[k];
        label idx = -1;
        forAll(availableNames, i)
        {
            if (availableNames[i] == name)
            {
                idx = i;
                break;
            }
        }

        if (idx < 0)
        {
            Info<< " " << name << "=<unknown>";
        }
        else
        {
            Info<< " " << name << "=" << values[idx];
        }
    }
    Info<< nl;
}

} // End namespace Foam
