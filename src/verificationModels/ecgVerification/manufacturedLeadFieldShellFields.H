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

Description
    Conductivity fields of the sphericalShell reference of
    manufacturedLeadFieldECGVerifier: 0/ConductivityIntracellular and
    0/ConductivityExtracellular (Gi, Ge of the transmission solution) and
    0/bodyAndOrgansConductivity (its sigmaB; the outer-bath formula in the
    myocardium, where the conductor does not read it), at the cell centres
    of the undecomposed heart+torso mesh. Each cell takes the formula of its
    cellZone (manufacturedLeadFieldShellZones).

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedLeadFieldShellFields_H
#define manufacturedLeadFieldShellFields_H

#include "volFields.H"
#include "bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H"

namespace Foam
{

//- Zone of the transmission solution for every cell, from the cellZones
//  bloodPool, myocardium and bath of the mesh (fatal if a cell is in none).
//  A cell centre may lie outside the radial range of its zone, where a
//  polyhedral interface cuts the sphere, by at most the cell size cbrt(V)
//  (the zone formulas extend smoothly across the spheres).
inline labelList manufacturedLeadFieldShellZones(const fvMesh& mesh)
{
    using namespace shellTransmissionMMS;

    const Params p = defaultParams();
    const wordList names({"bloodPool", "myocardium", "bath"});
    const label zoneIds[3] = {bloodPool, myocardium, bath};

    labelList zone(mesh.nCells(), -1);
    forAll(names, i)
    {
        const label zoneI = mesh.cellZones().findZoneID(names[i]);
        if (zoneI < 0)
        {
            FatalErrorInFunction
                << "The sphericalShell solution needs the cellZone "
                << names[i] << "." << exit(FatalError);
        }
        for (const label cellI : mesh.cellZones()[zoneI])
        {
            zone[cellI] = zoneIds[i];
        }
    }

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V().field();
    label unassigned = 0;
    label misplaced = 0;
    forAll(zone, cellI)
    {
        if (zone[cellI] < 0)
        {
            ++unassigned;
            continue;
        }
        const scalar r = mag(C[cellI]);
        const scalar rMin = zone[cellI] == bloodPool ? 0 : (zone[cellI] == myocardium ? p.r1 : p.r2);
        const scalar rMax = zone[cellI] == bloodPool ? p.r1 : (zone[cellI] == myocardium ? p.r2 : GREAT);
        const scalar size = Foam::cbrt(V[cellI]);
        if (r < rMin - size || r > rMax + size)
        {
            ++misplaced;
        }
    }
    reduce(unassigned, sumOp<label>());
    reduce(misplaced, sumOp<label>());
    if (unassigned > 0 || misplaced > 0)
    {
        FatalErrorInFunction
            << unassigned << " cells lie in no sphericalShell cellZone and "
            << misplaced << " cells lie farther than their own size outside "
            << "the radial range of their zone." << exit(FatalError);
    }

    return zone;
}


inline void writeManufacturedLeadFieldShellFields(const fvMesh& mesh)
{
    using namespace shellTransmissionMMS;

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "The sphericalShell lead-field conductivities are written on "
            << "the undecomposed case; run in serial." << exit(FatalError);
    }

    const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);
    const Params p = defaultParams();
    const vectorField& C = mesh.C().primitiveField();
    const unsigned int oldPrecision = IOstream::defaultPrecision(17);

    const auto io = [&mesh](const word& name)
    {
        return IOobject(name, "0", mesh, IOobject::NO_READ, IOobject::NO_WRITE);
    };

    const auto copyInternalToPatches = [&mesh](auto& field)
    {
        forAll(mesh.boundary(), patchI)
        {
            if (polyPatch::constraintType(mesh.boundaryMesh()[patchI].type()))
            {
                continue;
            }
            field.boundaryFieldRef()[patchI] ==
                field.boundaryField()[patchI].patchInternalField();
        }
    };

    for (label i = 0; i < 2; ++i)
    {
        volSymmTensorField G
        (
            io(i == 0 ? "ConductivityIntracellular" : "ConductivityExtracellular"),
            mesh,
            dimensionedSymmTensor(conductivityDim, Zero),
            "calculated"
        );
        forAll(C, cellI)
        {
            G[cellI] = i == 0 ? Gi(C[cellI], p) : Ge(C[cellI], p);
        }
        copyInternalToPatches(G);
        G.write();
        Info<< "Wrote " << G.name() << " in 0." << endl;
    }

    volScalarField sigma
    (
        io("bodyAndOrgansConductivity"),
        mesh,
        dimensionedScalar(conductivityDim, Zero),
        "calculated"
    );
    const labelList zone(manufacturedLeadFieldShellZones(mesh));
    forAll(C, cellI)
    {
        sigma[cellI] =
            sigmaB(C[cellI], zone[cellI] == myocardium ? label(bath) : zone[cellI], p);
    }
    copyInternalToPatches(sigma);
    sigma.write();
    Info<< "Wrote " << sigma.name() << " in 0." << endl;

    IOstream::defaultPrecision(oldPrecision);
}

} // End namespace Foam

#endif

// ************************************************************************* //
