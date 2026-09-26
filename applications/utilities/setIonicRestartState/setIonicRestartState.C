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

Application
    setIonicRestartState

Description
    Writes the per-cell ionic restart state file (<Model>State) that the
    myocardium reads at start-up. Each cell receives the state of its
    ionicHeterogeneity namedRegions region, read from one single-cell
    restart state file per region. Cells inside a blend transition receive
    the same weighted mix that is used for the region constants.

Usage
    setIonicRestartState [-region name] [-dict file] [-parallel]

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "ionicHeterogeneity.H"
#include "restartStateIO.H"

using namespace Foam;

// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace
{

//- Myocardium coefficients: <myocardiumSolver>Coeffs
const dictionary& myocardiumCoeffs(const dictionary& electroDict)
{
    const word solverName(electroDict.get<word>("myocardiumSolver"));
    return electroDict.subDict(solverName + "Coeffs");
}


//- Expand a path and resolve it against the case root
fileName casePath(fileName name, const Time& runTime)
{
    name.expand();
    if (!name.isAbsolute())
    {
        name = runTime.globalPath()/name;
    }
    return name;
}


//- State vector of a single-cell restart state file
scalarField readSingleCellState
(
    const fileName& statePath,
    const word& model,
    label& stride
)
{
    std::ifstream is(statePath.c_str(), std::ios::binary);
    if (!is)
    {
        FatalErrorInFunction
            << "Cannot open state file " << statePath
            << exit(FatalError);
    }

    const restartStateIO::Header header =
        restartStateIO::readHeader(is, statePath);

    if (header.model != model)
    {
        FatalErrorInFunction
            << "State file " << statePath << " is for ionic model '"
            << header.model << "', but the case uses '" << model << "'."
            << exit(FatalError);
    }

    if (header.rows != 1)
    {
        FatalErrorInFunction
            << "State file " << statePath << " holds " << header.rows
            << " cells; a single-cell state file is required."
            << exit(FatalError);
    }

    if (stride >= 0 && header.stride != stride)
    {
        FatalErrorInFunction
            << "State file " << statePath << " has " << header.stride
            << " states per cell; the other region files have " << stride
            << "."
            << exit(FatalError);
    }
    stride = header.stride;

    scalarField state(stride);
    forAll(state, stateI)
    {
        state[stateI] = restartStateIO::readScalar(is, statePath);
        restartStateIO::checkValue(state[stateI], statePath);
    }

    return state;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * Main  * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write the per-cell ionic restart state (<Model>State) of the"
        " myocardium from one single-cell state file per ionicHeterogeneity"
        " namedRegions region."
    );
    argList::addOption
    (
        "dict",
        "file",
        "State map dictionary (default: system/setIonicRestartStateDict)"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Region definition: the case's own ionicHeterogeneity block

    const IOdictionary electroDict
    (
        IOobject
        (
            "electroProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& coeffs = myocardiumCoeffs(electroDict);
    const word model(coeffs.get<word>("ionicModel"));

    if (coeffs.found("cellZone"))
    {
        FatalErrorInFunction
            << "Myocardium cellZone subsets are not supported: the state"
            << " file must cover every cell of the mesh."
            << exit(FatalError);
    }

    const dictionary& hetDict = coeffs.subDict("ionicHeterogeneity");

    const word mode(hetDict.get<word>("mode"));
    if (mode != "namedRegions")
    {
        FatalErrorInFunction
            << "ionicHeterogeneity mode '" << mode << "' is not supported."
            << " Supported: namedRegions."
            << exit(FatalError);
    }

    const word fieldName(hetDict.get<word>("field"));
    const word transitionMode(hetDict.get<word>("transitionMode"));

    // Width and smoothing are read only for blend; hard does not use them
    scalar transitionWidth = 0;
    word smoothing = word::null;

    if (transitionMode == "blend")
    {
        transitionWidth = hetDict.get<scalar>("transitionWidth");
        smoothing = hetDict.get<word>("smoothing");
    }
    else if (transitionMode != "hard")
    {
        FatalErrorInFunction
            << "ionicHeterogeneity transitionMode '" << transitionMode
            << "' is not supported. Supported: blend, hard."
            << exit(FatalError);
    }

    const List<ionicHeterogeneity::NamedFieldRegion> regions =
        ionicHeterogeneity::parseNamedFieldRegions(hetDict.subDict("regions"));

    wordList regionNames(regions.size());
    forAll(regions, regionI)
    {
        regionNames[regionI] = regions[regionI].name;
    }

    // State map: one single-cell state file per region

    const fileName dictPath
    (
        casePath
        (
            args.getOrDefault<fileName>
            (
                "dict",
                fileName("system/setIonicRestartStateDict")
            ),
            runTime
        )
    );

    IFstream mapStream(dictPath);
    if (!mapStream.good())
    {
        FatalErrorInFunction
            << "Cannot read state map dictionary " << dictPath
            << exit(FatalError);
    }
    const dictionary mapDict(mapStream);
    const dictionary& regionStates = mapDict.subDict("regionStates");

    for (const word& key : regionStates.toc())
    {
        if (!regionNames.found(key))
        {
            FatalErrorInFunction
                << "regionStates entry '" << key << "' is not an"
                << " ionicHeterogeneity region. Regions: " << regionNames
                << exit(FatalError);
        }
    }

    Info<< "Ionic model: " << model << nl
        << "Region field: " << fieldName << nl
        << "Region states:" << nl;

    label stride = -1;
    PtrList<scalarField> regionState(regions.size());

    forAll(regions, regionI)
    {
        if (!regionStates.found(regionNames[regionI]))
        {
            FatalErrorInFunction
                << "No state file for region '" << regionNames[regionI]
                << "' in " << dictPath << "."
                << exit(FatalError);
        }

        const fileName statePath
        (
            casePath
            (
                regionStates.get<fileName>(regionNames[regionI]),
                runTime
            )
        );

        regionState.set
        (
            regionI,
            new scalarField(readSingleCellState(statePath, model, stride))
        );

        Info<< "    " << regionNames[regionI] << ": " << statePath << nl;
    }

    // Per-cell state: region weights applied to the region states

    const volScalarField field
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    PtrList<scalarField> states(mesh.nCells());
    labelList dominantCount(regions.size(), 0);

    forAll(field, cellI)
    {
        const scalar rawT = field[cellI];

        if (rawT < -SMALL || rawT > 1 + SMALL)
        {
            FatalErrorInFunction
                << "Field '" << fieldName << "' value " << rawT
                << " at cell " << cellI << " is outside [0, 1]."
                << exit(FatalError);
        }

        const scalar t = min(max(rawT, scalar(0)), scalar(1));

        const List<ionicHeterogeneity::NamedRegionWeight> weights =
            ionicHeterogeneity::namedRegionWeightsAt
            (
                t, regions, transitionWidth, smoothing, transitionMode
            );

        scalarField cellState(stride, 0.0);
        label dominantI = -1;
        scalar dominantWeight = -1;

        forAll(weights, wI)
        {
            const label regionI = regionNames.find(weights[wI].name);
            cellState += weights[wI].weight*regionState[regionI];

            if (weights[wI].weight > dominantWeight)
            {
                dominantWeight = weights[wI].weight;
                dominantI = regionI;
            }
        }

        ++dominantCount[dominantI];
        states.set(cellI, new scalarField(cellState));
    }

    // Write through the solver's own restart format

    const fileName statePath =
        restartStateIO::path(mesh, model + "State");
    mkDir(statePath.path());
    restartStateIO::writeStates(mesh, model, stride, states);

    Info<< nl << "Wrote " << model << "State at time " << runTime.timeName()
        << ": " << returnReduce(mesh.nCells(), sumOp<label>())
        << " cells, " << stride << " states per cell" << nl
        << "Cells per region (largest weight):" << nl;

    forAll(regions, regionI)
    {
        Info<< "    " << regionNames[regionI] << ": "
            << returnReduce(dominantCount[regionI], sumOp<label>()) << nl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
