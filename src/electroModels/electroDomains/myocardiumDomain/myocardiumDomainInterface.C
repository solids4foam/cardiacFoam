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

#include "myocardiumDomainInterface.H"
#include "myocardiumDomain.H"
#include "eikonalMyocardiumDomain.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "error.H"
#include "fvMeshSubset.H"
#include "volFields.H"

namespace Foam
{

namespace
{

word myocardiumSolverType(const dictionary& electroProperties)
{
    const word coeffDictName = electroProperties.dictName();

    return coeffDictName.endsWith("Coeffs")
        ? word(coeffDictName.substr(0, coeffDictName.size() - 6))
        : coeffDictName;
}


scalarField readTransmuralDistance
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    const dictionary& heterogeneityDict
)
{
    scalarField fullValues;

    const word mode = heterogeneityDict.lookupOrDefault<word>("mode", "transmuralBands");

    if (mode == "cellZoneRegions")
    {
        fullValues.setSize(mesh.nCells(), -1.0);
        const dictionary& regionsDict = heterogeneityDict.subDict("regions");

        label regionIndex = 0;
        forAllConstIter(dictionary, regionsDict, iter)
        {
            if (iter().isDict())
            {
                const dictionary& regionDict = iter().dict();
                if (regionDict.found("cellZone"))
                {
                    const word zoneName(regionDict.lookup("cellZone"));
                    const label zoneId = mesh.cellZones().findZoneID(zoneName);

                    if (zoneId < 0)
                    {
                        FatalErrorInFunction
                            << "ionicHeterogeneity region '" << iter().keyword()
                            << "' specifies cellZone '" << zoneName
                            << "' but it does not exist on mesh '"
                            << mesh.name() << "'."
                            << exit(FatalError);
                    }

                    const labelList& zoneCells = mesh.cellZones()[zoneId];
                    forAll(zoneCells, i)
                    {
                        fullValues[zoneCells[i]] = scalar(regionIndex);
                    }
                }
                regionIndex++;
            }
        }
    }
    else
    {
        const word fieldName =
            heterogeneityDict.lookupOrDefault<word>("field", "t");

        // Construct temporary field to read values
        const volScalarField transmuralField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        fullValues = transmuralField.primitiveField();
    }

    if (!electroProperties.found("cellZone"))
    {
        return fullValues;
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = mesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << mesh.name() << "'."
            << exit(FatalError);
    }

    fvMeshSubset subset(mesh);
    subset.setCellSubset(mesh.cellZones()[zoneId]);

    const labelUList& cellMap = subset.cellMap();
    scalarField mappedValues(cellMap.size(), 0.0);

    forAll(cellMap, subCellI)
    {
        mappedValues[subCellI] = fullValues[cellMap[subCellI]];
    }

    return mappedValues;
}

} // End anonymous namespace


autoPtr<myocardiumDomainInterface> myocardiumDomainInterface::New
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    autoPtr<ionicModel>& ionicModelPtr,
    autoPtr<electroVerificationModel>& verificationModelPtr,
    scalar initialDeltaT
)
{
    const word solverType = myocardiumSolverType(electroProperties);

    if (solverType == "eikonalSolver")
    {
        ionicModelPtr.clear();
        verificationModelPtr.clear();

        return autoPtr<myocardiumDomainInterface>
        (
            new eikonalMyocardiumDomain(mesh, electroProperties)
        );
    }

    ionicModelPtr =
        ionicModel::New
        (
            electroProperties,
            myocardiumDomain::configuredCellCount(mesh, electroProperties),
            initialDeltaT
        );

    if (electroProperties.found("ionicHeterogeneity"))
    {
        const dictionary& heterogeneityDict =
            electroProperties.subDict("ionicHeterogeneity");

        if
        (
            heterogeneityDict.found("mode")
         || heterogeneityDict.found("field")
        )
        {
            const scalarField transmuralDistance =
                readTransmuralDistance(mesh, electroProperties, heterogeneityDict);

            ionicModelPtr->configureIonicHeterogeneity
            (
                transmuralDistance,
                heterogeneityDict
            );
        }

        if (heterogeneityDict.found("apexBaseBands"))
        {
            const dictionary& abDict =
                heterogeneityDict.subDict("apexBaseBands");

            const scalarField longitudinalDist =
                readTransmuralDistance(mesh, electroProperties, abDict);

            ionicModelPtr->configureApexBaseBandsHeterogeneity
            (
                longitudinalDist,
                abDict
            );
        }
    }

    verificationModelPtr =
        electroVerificationModel::New(electroProperties);

    return autoPtr<myocardiumDomainInterface>
    (
        myocardiumDomain::New
        (
            mesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModelPtr(),
            verificationModelPtr.get()
        ).ptr()
    );
}

} // End namespace Foam

// ************************************************************************* //
