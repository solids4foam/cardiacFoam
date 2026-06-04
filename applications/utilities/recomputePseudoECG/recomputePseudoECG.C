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
    recomputePseudoECG

Description
    Recompute the pseudo-ECG from stored Vm and conductivity fields without
    re-running the full simulation. Reads electrode positions from
    constant/electroProperties, then iterates over every stored time
    directory, applies the Gima-Rudy dipole integral:

        phi(P) = sum_c  [ (sigma_c . grad(Vm)_c) . (x_c - P) * V_c
                          / |x_c - P|^3 ]

    and writes postProcessing/pseudoECG.dat (or a user-specified file).

    Use this utility after changing electrode positions in electroProperties
    so you do not need to re-run the 3-D simulation.

Usage
    recomputeECG

    Options:
      -output   <path>   Override output file path
                         (default: postProcessing/pseudoECG.dat)
      -vmField  <name>   Name of the voltage field to read (default: Vm)
      -sigmaField <name> Name of the conductivity tensor field
                         (default: Diffusivity)
      -tStart   <scalar> Skip time directories earlier than this value
      -tEnd     <scalar> Skip time directories later than this value

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"

using namespace Foam;

List<word> readElectrodeNames
(
    const dictionary& electrodeDict
)
{
    List<word> names;
    forAllConstIter(dictionary, electrodeDict, iter)
    {
        names.append(iter().keyword());
    }
    return names;
}


List<vector> readElectrodePositions
(
    const dictionary& electrodeDict,
    const List<word>& names
)
{
    List<vector> positions(names.size());
    forAll(names, i)
    {
        positions[i] = electrodeDict.get<vector>(names[i]);
    }
    return positions;
}


const dictionary& findECGDict
(
    const dictionary& coeffs
)
{
    if (!coeffs.found("ecgDomains"))
    {
        FatalErrorInFunction
            << "Cannot find 'ecgDomains' in active myocardium coefficients "
            << " of constant/electroProperties."
            << exit(FatalError);
    }

    const dictionary& ecgDomains = coeffs.subDict("ecgDomains");
    const dictionary* fallbackDomainPtr = nullptr;

    forAllConstIter(dictionary, ecgDomains, iter)
    {
        if (iter().isDict() && iter().keyword() != "electrodePositions")
        {
            const dictionary& ecgDict = iter().dict();
            if (!fallbackDomainPtr)
            {
                fallbackDomainPtr = &ecgDict;
            }

            if (ecgDict.found("electrodePositions"))
            {
                return ecgDict;
            }
        }
    }

    if (fallbackDomainPtr && ecgDomains.found("electrodePositions"))
    {
        return *fallbackDomainPtr;
    }

    FatalErrorInFunction
        << "No ECG domain with 'electrodePositions' found inside ecgDomains."
        << exit(FatalError);

    return ecgDomains;
}


const dictionary& activeMyocardiumCoeffs
(
    const IOdictionary& electroProperties,
    word& solverType
)
{
    solverType =
        electroProperties.lookupOrDefault<word>
        (
            "myocardiumSolver",
            "monodomainSolver"
        );

    const word coeffsName(solverType + "Coeffs");

    if (electroProperties.found(coeffsName))
    {
        return electroProperties.subDict(coeffsName);
    }

    return electroProperties;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Recompute pseudo-ECG from stored Vm and conductivity fields.\n"
        "Electrode positions are read from constant/electroProperties.\n"
        "Output: postProcessing/pseudoECG.dat"
    );

    timeSelector::addOptions();

    argList::addOption
    (
        "output",
        "fileName",
        "Output file path (default: postProcessing/pseudoECG.dat)"
    );
    argList::addOption
    (
        "vmField",
        "word",
        "Vm field name to read from each time directory (default: Vm)"
    );
    argList::addOption
    (
        "sigmaField",
        "word",
        "Conductivity tensor field name in 0/ "
        "(default: Diffusivity)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    if (timeDirs.empty())
    {
        FatalErrorInFunction
            << "No time directories selected. Run from the case root "
            << "or use -time to specify a range."
            << exit(FatalError);
    }

    #include "createMesh.H"

    const word vmFieldName = args.getOrDefault<word>("vmField", "Vm");

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word myocardiumSolverType;
    const dictionary& myocardiumCoeffs =
        activeMyocardiumCoeffs(electroProperties, myocardiumSolverType);

    const word sigmaFieldName =
        args.getOrDefault<word>
        (
            "sigmaField",
            word("Diffusivity")
        );

    const dictionary& ecgDict = findECGDict(myocardiumCoeffs);

    const dictionary& ecgDomains = myocardiumCoeffs.subDict("ecgDomains");
    const dictionary* electrodeDictPtr = ecgDict.findDict("electrodePositions");

    if (!electrodeDictPtr)
    {
        electrodeDictPtr = ecgDomains.findDict("electrodePositions");
    }

    if (!electrodeDictPtr)
    {
        FatalErrorInFunction
            << "No electrodePositions dictionary found for the selected "
            << "pseudo-ECG domain."
            << exit(FatalError);
    }

    const dictionary& electrodeDict = *electrodeDictPtr;
    const List<word>   electrodeNames     = readElectrodeNames(electrodeDict);
    const List<vector> electrodePositions =
        readElectrodePositions(electrodeDict, electrodeNames);
    const label nElectrodes = electrodeNames.size();

    const scalar sigmaE =
        ecgDict.lookupOrDefault<scalar>("sigmaExtracellular", 0.0);
    const scalar normFactor =
        (sigmaE > VSMALL)
      ? 1.0 / (4.0 * constant::mathematical::pi * sigmaE)
      : 1.0;

    Info<< "Electrode positions:" << nl;
    forAll(electrodeNames, i)
    {
        Info<< "  " << electrodeNames[i] << " = " << electrodePositions[i] << nl;
    }
    Info<< endl;

    runTime.setTime(timeDirs[0], 0);
    mesh.readUpdate();

    volTensorField conductivity
    (
        IOobject
        (
            sigmaFieldName,
            "0",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "zero",
            pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume),
            tensor::zero
        )
    );

    if (!conductivity.headerOk())
    {
        Info<< "Conductivity tensor field '" << sigmaFieldName
            << "' not found in 0/, using value from active "
            << myocardiumSolverType << " coefficients." << nl << endl;

        conductivity = dimensionedTensor
        (
            dimensionedSymmTensor
            (
                sigmaFieldName,
                pow3(dimTime) * sqr(dimCurrent)/(dimMass*dimVolume),
                myocardiumCoeffs
            ) & tensor(I)
        );
    }

    const tensorField& sigma = conductivity.primitiveField();
    const scalarField& cellVolumes  = mesh.V();
    const vectorField& cellCentres  = mesh.C().primitiveField();

    Info<< "Read conductivity tensor field '" << sigmaFieldName
        << "' with " << sigma.size() << " cells." << nl << endl;

    const fileName outputPath
    (
        args.getOrDefault<fileName>
        (
            "output",
            runTime.globalPath()/"postProcessing"/"pseudoECG.dat"
        )
    );

    if (Pstream::master())
    {
        mkDir(outputPath.path());
    }

    autoPtr<OFstream> osPtr;
    if (Pstream::master())
    {
        osPtr.reset(new OFstream(outputPath));
        OFstream& os = osPtr();
        os << "# time";
        forAll(electrodeNames, i)
        {
            os << "  " << electrodeNames[i];
        }
        os << nl;
    }

    Info<< "Processing " << timeDirs.size() << " time steps..." << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();

        const scalar t = runTime.value();

        IOobject vmCheck
        (
            vmFieldName,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );
        if (!vmCheck.typeHeaderOk<volScalarField>(true))
        {
            Info<< "  Skipping t=" << t
                << " (no " << vmFieldName << " field)" << nl;
            continue;
        }

        volScalarField Vm
        (
            IOobject
            (
                vmFieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        const tmp<volVectorField> tgradVm = fvc::grad(Vm);
        const vectorField& gradVm = tgradVm().primitiveField();

        scalarField values(nElectrodes, 0.0);

        forAll(cellCentres, cellI)
        {
            const vector dipole =
                (sigma[cellI] & gradVm[cellI]) * cellVolumes[cellI];

            for (label eI = 0; eI < nElectrodes; ++eI)
            {
                const vector rVec = cellCentres[cellI] - electrodePositions[eI];
                const scalar r    = mag(rVec);

                if (r > VSMALL)
                {
                    values[eI] += (dipole & rVec) / (r * r * r);
                }
            }
        }

        for (label eI = 0; eI < nElectrodes; ++eI)
        {
            reduce(values[eI], sumOp<scalar>());
        }

        forAll(values, eI)
        {
            values[eI] *= normFactor;
        }

        if (Pstream::master())
        {
            OFstream& os = osPtr();
            os << t;
            forAll(values, eI)
            {
                os << "  " << values[eI];
            }
            os << nl;
        }

        if (timeI % 50 == 0 || timeI == timeDirs.size() - 1)
        {
            Info<< "  t = " << t
                << "  (" << timeI + 1 << "/" << timeDirs.size() << ")" << nl;
        }
    }

    Info<< nl << "Written: " << outputPath << nl;
    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
