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

#include "torsoCoupling.H"
#include "dictionaryKeyCheck.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::torsoCoupling::type Foam::torsoCoupling::read
(
    const dictionary& electroProperties
)
{
    const bool hasTorso = electroProperties.isDict("torso");
    const bool hasCoupling = electroProperties.found("torsoCoupling");

    if (!hasTorso && !hasCoupling)
    {
        return type::none;
    }

    if (!hasCoupling)
    {
        FatalIOErrorInFunction(electroProperties)
            << "electroProperties has a torso block but no torsoCoupling. "
            << "Set torsoCoupling bath (phiE of the heart and the torso "
            << "solved every step, bidomainSolver only) or torsoCoupling "
            << "leadField (lead fields precomputed by buildLeadFields, "
            << "monodomainSolver or bidomainSolver)."
            << exit(FatalIOError);
    }

    const word value
    (
        dictionaryKeyCheck::getWord
        (
            electroProperties,
            "torsoCoupling",
            word::null,
            wordList({"bath", "leadField"})
        )
    );

    if (!hasTorso)
    {
        FatalIOErrorInFunction(electroProperties)
            << "electroProperties sets torsoCoupling " << value << " but has "
            << "no torso block. Add" << nl
            << "torso" << nl
            << "{" << nl
            << "    cellZones         (...);" << nl
            << "    conductivityField <name>;" << nl
            << "    interfaceConductivityInterpolation "
            << "distanceWeightedHarmonic;" << nl
            << "    ..." << nl
            << "}" << exit(FatalIOError);
    }

    return value == "bath" ? type::bath : type::leadField;
}


Foam::torsoCoupling::type Foam::torsoCoupling::readFromCoeffs
(
    const dictionary& myocardiumCoeffs
)
{
    const dictionary& electroProperties = myocardiumCoeffs.topDict();

    if (&electroProperties == &myocardiumCoeffs)
    {
        FatalIOErrorInFunction(myocardiumCoeffs)
            << "Dictionary " << myocardiumCoeffs.name() << " is not inside "
            << "electroProperties, which holds torso and torsoCoupling."
            << exit(FatalIOError);
    }

    return read(electroProperties);
}


Foam::word Foam::torsoCoupling::name(const type coupling)
{
    switch (coupling)
    {
        case type::bath: return "bath";
        case type::leadField: return "leadField";
        default: return "none";
    }
}


// ************************************************************************* //
