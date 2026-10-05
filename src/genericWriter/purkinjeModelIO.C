/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "purkinjeModelIO.H"
#include "OSspecific.H"
#include "Time.H"
#include "PstreamReduceOps.H"
#include <fstream>
#include <sstream>
#include "ionicVariableCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

purkinjeModelIO::ResolvedTokens purkinjeModelIO::filterTokens
(
    const wordList& userExport,
    const char* const stateNames[],
    int nStates,
    const char* const algNames[],
    int nAlg
)
{
    ResolvedTokens res;
    DynamicList<word> networkTokens(userExport.size());
    DynamicList<word> ionicTokens(userExport.size());
    DynamicList<word> unknownTokens;

    for (const word& var : userExport)
    {
        if (var == "Vm" || var == "Iion" || var == "activationTime" ||
            var == "IcouplingSource" || var == "IcouplingCurrent")
        {
            networkTokens.append(var);
        }
        else
        {
            bool isVmDummy = false;
            label stateIdx = -1;
            label algIdx = -1;
            label rateIdx = -1;

            bool resolved = ionicVariableCompatibility::resolveVariable
            (
                var,
                stateNames,
                nStates,
                algNames,
                nAlg,
                isVmDummy,
                stateIdx,
                algIdx,
                rateIdx
            );

            if (resolved && (stateIdx >= 0 || algIdx >= 0))
            {
                ionicTokens.append(var);
            }
            else
            {
                unknownTokens.append(var);
            }
        }
    }

    res.networkTokens = networkTokens;
    res.ionicTokens = ionicTokens;
    res.unknownTokens = unknownTokens;
    return res;
}

autoPtr<OFstream> purkinjeModelIO::openTimeSeries
(
    const fileName& outDir,
    const word& filename,
    const wordList& columnNames,
    const scalar startTime
)
{
    mkDir(outDir);

    std::string header("# time");
    forAll(columnNames, i)
    {
        header += "  " + columnNames[i];
    }

    DynamicList<std::string> keptRows;
    if (startTime > 0 && isFile(outDir/filename))
    {
        // Rows are written with 8 significant digits.
        const scalar lastKept = startTime*(1 + 1e-6);

        std::ifstream is((outDir/filename).c_str());
        std::string line;
        if (std::getline(is, line) && line == header)
        {
            while (std::getline(is, line))
            {
                std::istringstream row(line);
                double t = 0;
                if ((row >> t) && t <= lastKept)
                {
                    keptRows.append(line);
                }
            }
        }
    }

    autoPtr<OFstream> osPtr(new OFstream(outDir/filename));

    OFstream& os = osPtr.ref();
    os.setf(std::ios::scientific);
    os.precision(8);

    os << header.c_str() << nl;
    forAll(keptRows, i)
    {
        os << keptRows[i].c_str() << nl;
    }

    return osPtr;
}


void purkinjeModelIO::readVTKSeries
(
    const fileName& seriesPath,
    const scalar startTime,
    scalarList& times,
    wordList& filenames
)
{
    DynamicList<scalar> keptTimes;
    DynamicList<word> keptFiles;

    std::ifstream is(seriesPath.c_str());
    std::string line;
    const std::string nameKey("\"name\": \"");
    const std::string timeKey("\"time\": ");
    while (std::getline(is, line))
    {
        const auto n = line.find(nameKey);
        const auto t = line.find(timeKey);
        if (n == std::string::npos || t == std::string::npos)
        {
            continue;
        }

        const auto nameStart = n + nameKey.size();
        const word name(line.substr(nameStart, line.find('"', nameStart) - nameStart));
        const scalar time = std::stod(line.substr(t + timeKey.size()));

        if (time <= startTime*(1 + 1e-6))
        {
            keptTimes.append(time);
            keptFiles.append(name);
        }
    }

    times = keptTimes;
    filenames = keptFiles;
}


void purkinjeModelIO::writeGlobalField(scalarGlobalIOField& field)
{
    bool writeGood = true;
    const Time& runTime = field.time();

    field.instance() = runTime.timeName();

    if (Pstream::master())
    {
        const fileName outputPath
        (
            runTime.globalPath()/field.instance()/field.name()
        );

        mkDir(outputPath.path());
        OFstream os
        (
            outputPath,
            IOstreamOption
            (
                runTime.writeFormat(),
                runTime.writeCompression()
            )
        );

        writeGood = os.good() && field.writeHeader(os) && field.writeData(os);

        if (writeGood)
        {
            IOobject::writeEndDivider(os);
        }
    }

    reduce(writeGood, andOp<bool>());

    if (!writeGood)
    {
        FatalErrorInFunction
            << "Failed writing " << field.name()
            << exit(FatalError);
    }
}


void purkinjeModelIO::writeRow
(
    OFstream& os,
    scalar time,
    const List<scalar>& values
)
{
    os << time;
    forAll(values, i)
    {
        os << "  " << values[i];
    }
    os << nl;
}


void purkinjeModelIO::writeVTKSeries
(
    const fileName&   seriesPath,
    const scalarList& times,
    const wordList&   filenames
)
{
    std::ofstream js(seriesPath);
    js.setf(std::ios::scientific);
    js.precision(8);

    js << "{\n"
       << "  \"file-series-version\": \"1.0\",\n"
       << "  \"files\": [\n";

    forAll(times, i)
    {
        js << "    { \"name\": \"" << filenames[i]
           << "\", \"time\": " << times[i] << " }";
        if (i < times.size() - 1) js << ",";
        js << "\n";
    }

    js << "  ]\n"
       << "}\n";
}


void purkinjeModelIO::writeVTK
(
    const fileName&    outDir,
    const word&        timeName,
    const label        timeStep,
    const pointField&  nodeLocations,
    const labelList&   edgeNodeA,
    const labelList&   edgeNodeB,
            const scalarField& Vm1D,
            const scalarField& Iion1D,
            const labelList&   pvjNodes,
            const scalarField& terminalSource,
            bool                includeVmIion,
            const PtrList<scalarField>& ionicFields,
            const wordList&             ionicFieldNames
        )
{
    if (nodeLocations.empty())
    {
        return;
    }

    mkDir(outDir);

    // Format filename with zero-padded timeStep index for correct ParaView sorting
    char buf[32];
    snprintf(buf, sizeof(buf), "purkinjeNetwork_%06d.vtk", int(timeStep));
    const fileName vtkPath(outDir / buf);

    OFstream os(vtkPath);
    os.setf(std::ios::scientific);
    os.precision(6);

    const label nNodes = nodeLocations.size();
    const label nEdges = edgeNodeA.size();

    // ---- VTK legacy ASCII header ----
    os  << "# vtk DataFile Version 2.0\n"
        << "Purkinje network t=" << timeName << "\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n";

    // ---- Node positions ----
    os  << "POINTS " << nNodes << " float\n";
    forAll(nodeLocations, i)
    {
        const point& p = nodeLocations[i];
        os  << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // ---- Edges as VTK line cells ("2 nodeA nodeB" per line) ----
    os  << "\nLINES " << nEdges << " " << 3*nEdges << "\n";
    forAll(edgeNodeA, eI)
    {
        os  << "2 " << edgeNodeA[eI] << " " << edgeNodeB[eI] << "\n";
    }

    // ---- Point data fields ----
    os  << "\nPOINT_DATA " << nNodes << "\n";

    if (includeVmIion)
    {
        // Vm natively in Volts
        os  << "SCALARS Vm_V float 1\n"
            << "LOOKUP_TABLE default\n";
        forAll(Vm1D, i)
        {
            os  << Vm1D[i] << "\n";
        }

        // Ionic current
        os  << "SCALARS Iion float 1\n"
            << "LOOKUP_TABLE default\n";
        forAll(Iion1D, i)
        {
            os  << Iion1D[i] << "\n";
        }
    }

    // Volumetric coupling source — non-zero only at PVJ nodes. Eikonal
    // couplers work on arrival time, not a volumetric source, so this is
    // meaningless (permanently zero) when there is no ionic model.
    if (includeVmIion && pvjNodes.size() && terminalSource.size() == pvjNodes.size())
    {
        scalarField pvjField(nNodes, 0.0);
        forAll(pvjNodes, k)
        {
            pvjField[pvjNodes[k]] = terminalSource[k];
        }

        os  << "SCALARS IcouplingSource_Am3 float 1\n"
            << "LOOKUP_TABLE default\n";
        forAll(pvjField, i)
        {
            os  << pvjField[i] << "\n";
        }
    }

    // ---- Dynamic Ionic Fields ----
    forAll(ionicFields, fI)
    {
        os  << "SCALARS " << ionicFieldNames[fI] << " float 1\n"
            << "LOOKUP_TABLE default\n";
        const scalarField& field = ionicFields[fI];
        forAll(field, i)
        {
            os << field[i] << "\n";
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
