/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    setTorsoOrganConductivityField

Description
    Create a scalar torso/bath conductivity field from named mesh cellZones.

    The utility reads system/setTorsoOrganConductivityFieldDict by default.
    Conductivities are assigned from the cellZones sub-dictionary. Unknown
    cellZone names are fatal. A default conductivity is optional and is only
    used for cells not covered by any configured zone.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

using namespace Foam;

namespace
{

const dimensionSet conductivityDim
(
    pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume)
);


void printAvailableCellZones(const fvMesh& mesh, Ostream& os)
{
    os << "Available cellZones:";

    if (mesh.cellZones().empty())
    {
        os << " none";
    }
    else
    {
        forAll(mesh.cellZones(), zoneI)
        {
            os << nl << "  " << mesh.cellZones()[zoneI].name()
               << " (" << mesh.cellZones()[zoneI].size() << " cells)";
        }
    }

    os << nl;
}


label countAssigned(const boolList& assigned)
{
    label n = 0;

    forAll(assigned, cellI)
    {
        if (assigned[cellI])
        {
            ++n;
        }
    }

    return n;
}

}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a torso/bath conductivity volScalarField from named cellZones."
    );
    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary path. Default: system/"
        "setTorsoOrganConductivityFieldDict"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const fileName dictPath =
    (
        args.found("dict")
      ? args.get<fileName>("dict")
      : runTime.system()/"setTorsoOrganConductivityFieldDict"
    );

    Info<< nl
        << "========== setTorsoOrganConductivityField =========="
        << nl << "Reading " << dictPath << nl << endl;

    IOdictionary dict
    (
        IOobject
        (
            dictPath.name(),
            dictPath.path(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word fieldName
    (
        dict.lookupOrDefault<word>("fieldName", "bodyAndOrgansConductivity")
    );

    if (!dict.found("cellZones"))
    {
        FatalErrorInFunction
            << "Missing required 'cellZones' sub-dictionary in "
            << dictPath << "." << nl
            << "Each entry must map a mesh cellZone name to a scalar "
            << "conductivity value."
            << exit(FatalError);
    }

    const dictionary& zoneDict = dict.subDict("cellZones");
    const wordList zoneNames(zoneDict.toc());

    if (zoneNames.empty())
    {
        FatalErrorInFunction
            << "'cellZones' in " << dictPath
            << " is empty; no conductivity mapping was configured."
            << exit(FatalError);
    }

    const bool haveDefault = dict.found("defaultSigma");
    const scalar defaultSigma =
        haveDefault ? dict.get<scalar>("defaultSigma") : 0.0;

    volScalarField conductivity
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", conductivityDim, 0.0),
        "zeroGradient"
    );

    boolList assigned(mesh.nCells(), false);

    Info<< "Writing field '" << fieldName << "'" << nl << endl;
    printAvailableCellZones(mesh, Info);
    Info<< nl << "Configured conductivity mapping:" << nl;

    forAll(zoneNames, zoneI)
    {
        const word& zoneName = zoneNames[zoneI];
        const label zoneId = mesh.cellZones().findZoneID(zoneName);

        if (zoneId < 0)
        {
            FatalErrorInFunction
                << "Configured cellZone '" << zoneName
                << "' was not found on mesh '" << mesh.name() << "'."
                << nl
                << "There is no valid mapping for this key."
                << nl << nl
                << "Check the cellZone names in "
                << "constant/polyMesh/cellZones or in the converted mesh."
                << nl
                << exit(FatalError);
        }

        const scalar sigma = zoneDict.get<scalar>(zoneName);
        const cellZone& zone = mesh.cellZones()[zoneId];

        label assignedHere = 0;

        forAll(zone, zoneCellI)
        {
            const label cellI = zone[zoneCellI];

            if (assigned[cellI])
            {
                FatalErrorInFunction
                    << "Cell " << cellI << " is assigned by more than one "
                    << "configured organ cellZone. Overlapping conductivity "
                    << "zones are ambiguous."
                    << exit(FatalError);
            }

            conductivity[cellI] = sigma;
            assigned[cellI] = true;
            ++assignedHere;
        }

        Info<< "  " << zoneName << ": sigma=" << sigma
            << ", cells=" << assignedHere << nl;
    }

    const label nAssignedBeforeDefault = countAssigned(assigned);
    const label nUnassigned = mesh.nCells() - nAssignedBeforeDefault;

    if (nUnassigned > 0)
    {
        if (!haveDefault)
        {
            FatalErrorInFunction
                << nUnassigned << " of " << mesh.nCells()
                << " cells were not covered by configured cellZones and no "
                << "'defaultSigma' was provided." << nl
                << "Add the missing organ cellZone mapping or explicitly set "
                << "defaultSigma if fallback assignment is intended."
                << exit(FatalError);
        }

        forAll(assigned, cellI)
        {
            if (!assigned[cellI])
            {
                conductivity[cellI] = defaultSigma;
                assigned[cellI] = true;
            }
        }

        Info<< "  defaultSigma: sigma=" << defaultSigma
            << ", cells=" << nUnassigned
            << " (fallback for cells not in configured zones)" << nl;
    }
    else if (haveDefault)
    {
        Info<< "  defaultSigma: configured but not used; all cells were "
            << "covered by named cellZones." << nl;
    }

    conductivity.correctBoundaryConditions();

    Info<< nl
        << "Assigned cells: " << countAssigned(assigned)
        << " / " << mesh.nCells() << nl
        << "min(" << fieldName << ") = " << gMin(conductivity)
        << ", max(" << fieldName << ") = " << gMax(conductivity)
        << nl << endl;

    conductivity.write();

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
