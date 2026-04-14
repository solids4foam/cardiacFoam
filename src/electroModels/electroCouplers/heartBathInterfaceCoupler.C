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

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  Stub for heart-bath interface coupler.
\*---------------------------------------------------------------------------*/

#include "heartBathInterfaceCoupler.H"

#include "DynamicList.H"
#include "UIndirectList.H"
#include "addToRunTimeSelectionTable.H"
#include "electroDomainInterface.H"
#include "fvPatch.H"

namespace Foam
{

namespace
{

scalar nearestFaceDistanceSqr
(
    const point& target,
    const vectorField& candidates,
    label& nearestFacei
)
{
    scalar bestDistSqr = GREAT;
    nearestFacei = -1;

    forAll(candidates, facei)
    {
        const scalar distSqr = magSqr(candidates[facei] - target);

        if (distSqr < bestDistSqr)
        {
            bestDistSqr = distSqr;
            nearestFacei = facei;
        }
    }

    return bestDistSqr;
}

bathCouplingEndpoint& requireBathCouplingDomain
(
    electroDomainInterface& secondaryDomain
)
{
    bathCouplingEndpoint* bathDomain =
        dynamic_cast<bathCouplingEndpoint*>(&secondaryDomain);

    if (!bathDomain)
    {
        FatalErrorInFunction
            << "Configured bath domain does not implement "
            << "bathCouplingEndpoint."
            << exit(FatalError);
    }

    return *bathDomain;
}

const electroStateProvider& requireHeartBathStateProvider
(
    tissueCouplingEndpoint& primaryDomain
)
{
    const electroStateProvider* heartStateProvider =
        dynamic_cast<const electroStateProvider*>(&primaryDomain);

    if (!heartStateProvider)
    {
        FatalErrorInFunction
            << "Configured heart-to-bath primary domain does not implement "
            << "electroStateProvider."
            << exit(FatalError);
    }

    return *heartStateProvider;
}

} // End anonymous namespace


defineTypeNameAndDebug(heartBathInterfaceCoupler, 0);
addToRunTimeSelectionTable
(
    ElectroDomainCoupler,
    heartBathInterfaceCoupler,
    dictionary
);


bathCouplingEndpoint& heartBathInterfaceCoupler::requireBathDomain
(
    electroDomainInterface& secondaryDomain
)
{
    return Foam::requireBathCouplingDomain(secondaryDomain);
}


const electroStateProvider& heartBathInterfaceCoupler::requireHeartStateProvider
(
    tissueCouplingEndpoint& primaryDomain
)
{
    return Foam::requireHeartBathStateProvider(primaryDomain);
}


heartBathInterfaceCoupler::heartBathInterfaceCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    ElectroDomainCoupler(primaryDomain, secondaryDomain),
    bathDomain_(requireBathDomain(secondaryDomain)),
    heartStateProvider_(requireHeartStateProvider(primaryDomain)),
    copyHeartPhiEToBath_
    (
        dict.lookupOrDefault<Switch>("copyHeartPhiEToBath", true)
    ),
    reportCouplingExtrema_
    (
        dict.lookupOrDefault<Switch>("reportCouplingExtrema", true)
    )
{
}


void heartBathInterfaceCoupler::preparePostPrimaryCoupling
(
    scalar t0,
    scalar dt
)
{
    (void)t0;
    (void)dt;

    bathDomain_.interfaceSourceField() =
        dimensionedScalar
        (
            "zero",
            bathDomain_.interfaceSourceField().dimensions(),
            0.0
        );
    bathDomain_.interfaceSourceField().correctBoundaryConditions();

    if (!copyHeartPhiEToBath_)
    {
        return;
    }

    const volScalarField* phiEPtr = heartStateProvider_.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "heartBathInterfaceCoupler requires the primary domain to "
            << "expose phiE. Select a bidomain primary solver or disable "
            << "copyHeartPhiEToBath in the coupling dictionary."
            << exit(FatalError);
    }

    volScalarField& bathPhiE = bathDomain_.bathPotentialField();
    const fvMesh& heartBaseMesh = heartStateProvider_.baseMesh();
    const scalarField& heartPhiE = phiEPtr->primitiveField();

    if
    (
        bathDomain_.subsetFaceMapPtr()
     && bathDomain_.subsetCellMapPtr()
     && heartStateProvider_.subsetCellMapPtr()
     && &bathDomain_.baseMesh() == &heartBaseMesh
    )
    {
        const label interfacePatchi = bathDomain_.interfacePatchIndex();

        if (interfacePatchi < 0)
        {
            FatalErrorInFunction
                << "Bath submesh coupling requires an interface patch."
                << exit(FatalError);
        }

        const labelUList& faceMap = *bathDomain_.subsetFaceMapPtr();
        const labelUList& bathCellMap = *bathDomain_.subsetCellMapPtr();
        const labelUList& heartCellMap = *heartStateProvider_.subsetCellMapPtr();
        const fvMesh& baseMesh = heartBaseMesh;
        const fvPatch& bathPatch = bathPhiE.mesh().boundary()[interfacePatchi];
        fvPatchScalarField& bathPatchField =
            bathPhiE.boundaryFieldRef()[interfacePatchi];
        scalarField mappedPatchValues(bathPatch.size(), 0.0);

        const label start = bathPatch.start();
        const labelUList& bathFaceCells = bathPatch.faceCells();
        const labelUList& owner = baseMesh.faceOwner();
        const labelUList& neighbour = baseMesh.faceNeighbour();
        labelList baseToHeartCell(baseMesh.nCells(), -1);

        forAll(heartCellMap, heartCellI)
        {
            baseToHeartCell[heartCellMap[heartCellI]] = heartCellI;
        }

        if (reportCouplingExtrema_)
        {
            Info<< "heartBathInterfaceCoupler: interface patch geometry size = "
                << bathPatch.size()
                << ", patch field size = " << bathPatchField.size()
                << nl;
        }

        forAll(mappedPatchValues, facei)
        {
            const label baseFace = faceMap[start + facei];

            if (baseFace >= baseMesh.nInternalFaces())
            {
                FatalErrorInFunction
                    << "Mapped bath interface face " << baseFace
                    << " is not an internal face on base mesh '"
                    << baseMesh.name() << "'."
                    << exit(FatalError);
            }

            const label baseBathCell = bathCellMap[bathFaceCells[facei]];
            const label ownerCell = owner[baseFace];
            const label neighbourCell = neighbour[baseFace];
            const label baseHeartCell =
                (ownerCell == baseBathCell ? neighbourCell : ownerCell);
            const label heartCell = baseToHeartCell[baseHeartCell];

            if (heartCell < 0)
            {
                FatalErrorInFunction
                    << "Base cell " << baseHeartCell
                    << " adjacent to the bath interface is not part of the "
                    << "heart submesh."
                    << exit(FatalError);
            }

            mappedPatchValues[facei] = heartPhiE[heartCell];
        }

        bathPatchField == mappedPatchValues;

        if (reportCouplingExtrema_)
        {
            Info<< "heartBathInterfaceCoupler: heart phiE min/max = "
                << gMin(heartPhiE) << " / " << gMax(heartPhiE)
                << ", mapped bath patch min/max = "
                << gMin(mappedPatchValues) << " / " << gMax(mappedPatchValues)
                << nl;
        }

        bathPhiE.correctBoundaryConditions();
        return;
    }

    if
    (
        !bathDomain_.subsetFaceMapPtr()
     && heartStateProvider_.subsetCellMapPtr()
     && bathDomain_.interfacePatchIndex() >= 0
     && &bathDomain_.baseMesh() == &heartBaseMesh
    )
    {
        const labelUList& heartCellMap = *heartStateProvider_.subsetCellMapPtr();
        const fvMesh& baseMesh = heartBaseMesh;
        const labelUList& owner = baseMesh.faceOwner();
        const labelUList& neighbour = baseMesh.faceNeighbour();
        const fvPatch& bathPatch = bathPhiE.mesh().boundary()[bathDomain_.interfacePatchIndex()];
        fvPatchScalarField& bathPatchField =
            bathPhiE.boundaryFieldRef()[bathDomain_.interfacePatchIndex()];
        scalarField mappedPatchValues(bathPatch.patch().size(), 0.0);
        labelList baseToHeartCell(baseMesh.nCells(), -1);
        DynamicList<label> heartInterfaceFaces;
        DynamicList<label> heartInterfaceCells;

        forAll(heartCellMap, heartCellI)
        {
            baseToHeartCell[heartCellMap[heartCellI]] = heartCellI;
        }

        for (label facei = 0; facei < baseMesh.nInternalFaces(); ++facei)
        {
            const bool ownerIsHeart = baseToHeartCell[owner[facei]] >= 0;
            const bool neighbourIsHeart = baseToHeartCell[neighbour[facei]] >= 0;

            if (ownerIsHeart == neighbourIsHeart)
            {
                continue;
            }

            heartInterfaceFaces.append(facei);
            heartInterfaceCells.append
            (
                ownerIsHeart
              ? baseToHeartCell[owner[facei]]
              : baseToHeartCell[neighbour[facei]]
            );
        }

        const vectorField baseInterfaceCentres
        (
            UIndirectList<vector>(baseMesh.faceCentres(), heartInterfaceFaces)
        );
        const vectorField& bathPatchCentres = bathPatch.Cf();

        if (reportCouplingExtrema_)
        {
            Info<< "heartBathInterfaceCoupler: real-region patch geometry size = "
                << bathPatch.patch().size()
                << ", patch field size = " << bathPatchField.size()
                << ", detected heart interface faces = "
                << heartInterfaceFaces.size() << nl;
        }

        forAll(mappedPatchValues, facei)
        {
            label nearestFacei = -1;
            const scalar distSqr = nearestFaceDistanceSqr
            (
                bathPatchCentres[facei],
                baseInterfaceCentres,
                nearestFacei
            );

            if (nearestFacei < 0 || distSqr > SMALL)
            {
                FatalErrorInFunction
                    << "Could not match bath interface face centre "
                    << bathPatchCentres[facei]
                    << " to a heart interface face on base mesh '"
                    << baseMesh.name() << "'. Best squared distance = "
                    << distSqr << "."
                    << exit(FatalError);
            }

            mappedPatchValues[facei] = heartPhiE[heartInterfaceCells[nearestFacei]];
        }

        bathPatchField == mappedPatchValues;

        if (reportCouplingExtrema_)
        {
            Info<< "heartBathInterfaceCoupler: heart phiE min/max = "
                << gMin(heartPhiE) << " / " << gMax(heartPhiE)
                << ", mapped bath patch min/max = "
                << gMin(mappedPatchValues) << " / " << gMax(mappedPatchValues)
                << nl;
        }

        bathPhiE.correctBoundaryConditions();
        return;
    }

    if (&phiEPtr->mesh() != &bathPhiE.mesh())
    {
        return;
    }

    bathPhiE.primitiveFieldRef() = phiEPtr->primitiveField();
    bathPhiE.correctBoundaryConditions();

    if (reportCouplingExtrema_)
    {
        Info<< "heartBathInterfaceCoupler: heart phiE min/max = "
            << gMin(heartPhiE) << " / " << gMax(heartPhiE)
            << ", copied bath phiE min/max = "
            << gMin(bathPhiE.primitiveField()) << " / "
            << gMax(bathPhiE.primitiveField()) << nl;
    }
}

} // End namespace Foam

// ************************************************************************* //
