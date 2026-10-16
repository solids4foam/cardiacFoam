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

#include "myocardiumSolver.H"
#include "error.H"
#include "nonOrthogonalCorrectorLoop.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(myocardiumSolver, 0);
defineRunTimeSelectionTable(myocardiumSolver, dictionary);


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

autoPtr<myocardiumSolver> myocardiumSolver::New
(
    const fvMesh& mesh,
    const fvMesh& supportMesh,
    const fvMeshSubset* meshSubsetPtr,
    const word& solverType,
    const dictionary& coeffs
)
{


    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            coeffs,
            "myocardiumSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<myocardiumSolver>
    (
        ctorPtr(mesh, supportMesh, meshSubsetPtr, coeffs)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myocardiumSolver::bindExternalPhiE
(
    const volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    (void)phiE;
    (void)heartCellMap;

    FatalErrorInFunction
        << "The selected myocardiumSolver cannot bind an external phiE field."
        << exit(FatalError);
}


void myocardiumSolver::unbindExternalPhiE()
{
    FatalErrorInFunction
        << "The selected myocardiumSolver cannot unbind an external phiE field."
        << exit(FatalError);
}


void myocardiumSolver::solveDiffusionImplicit
(
    electroVolumeFieldDomain& domain,
    scalar dt,
    pimpleControl& pimple
)
{
    while (pimple.loop())
    {
        correctNonOrthogonalLoop
        (
            pimple,
            [&]() { solveDiffusionImplicit(domain, dt); }
        );
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace
{

//- Gather the support-mesh file fieldName through cellMap; false if there is none
bool gatherMappedCellValues
(
    const word& fieldName,
    const fvMesh& supportMesh,
    const labelUList& cellMap,
    scalarField& subValues
)
{
    IOobject fullIO
    (
        fieldName,
        supportMesh.time().timeName(),
        supportMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!fullIO.typeHeaderOk<volScalarField>(true))
    {
        return false;
    }

    const volScalarField fullField(fullIO, supportMesh);
    const scalarField& fullValues = fullField.primitiveField();

    subValues.setSize(cellMap.size());
    forAll(cellMap, subCelli)
    {
        subValues[subCelli] = fullValues[cellMap[subCelli]];
    }

    return true;
}

} // End anonymous namespace


IOobject::readOption myocardiumFieldReadOption
(
    const fvMeshSubset* meshSubsetPtr
)
{
    return
    (
        meshSubsetPtr && meshSubsetPtr->hasSubMesh()
      ? IOobject::NO_READ
      : IOobject::READ_IF_PRESENT
    );
}


bool readMappedCellField
(
    volScalarField& subField,
    const fvMesh& supportMesh,
    const labelUList& cellMap
)
{
    scalarField values;

    if (!gatherMappedCellValues(subField.name(), supportMesh, cellMap, values))
    {
        return false;
    }

    subField.primitiveFieldRef() = values;
    subField.correctBoundaryConditions();

    if
    (
        gatherMappedCellValues
        (
            subField.name() + "_0", supportMesh, cellMap, values
        )
    )
    {
        volScalarField& subField0 = subField.oldTime();
        subField0.primitiveFieldRef() = values;
        subField0.correctBoundaryConditions();
        subField0.oldTime();
    }

    return true;
}

} // End namespace Foam

// ************************************************************************* //
