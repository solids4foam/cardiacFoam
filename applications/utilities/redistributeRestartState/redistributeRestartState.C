/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 cardiacFoam authors
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

Application
    redistributeRestartState

Description
    Moves the restart state files (<model>State) that cardiacFoam's ionic and
    active-tension models write across the processor boundary, which
    decomposePar and reconstructPar leave alone.

    By default it reconstructs: for every selected time, the files found in
    processor0/<time>/<region>/ are gathered with each processor's
    cellProcAddressing into <time>/<region>/<model>State of the undecomposed
    case. With -decompose it scatters the undecomposed files into the
    processor directories instead.

    A state file holds one row per cell. Files whose row count differs from
    the processor's cell count are reported and skipped.

Usage
    redistributeRestartState [-decompose] [-region <name>] [time options]
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "polyMesh.H"
#include "labelIOList.H"
#include "OSspecific.H"
#include "fileOperation.H"
#include "restartStateIO.H"

#include <fstream>

using namespace Foam;

namespace
{

//- The <model>State files in a directory, sorted
wordList stateFiles(const fileName& dir)
{
    wordList found;
    if (!isDir(dir))
    {
        return found;
    }
    for (const fileName& f : readDir(dir, fileName::FILE))
    {
        if (f.ends_with("State"))
        {
            found.append(word(f));
        }
    }
    Foam::sort(found);
    return found;
}


//- Time directory of a region, in the case or a processor directory
fileName regionTimeDir
(
    const fileName& casePath,
    const word& timeName,
    const word& regionDir
)
{
    return regionDir.empty() ? casePath/timeName : casePath/timeName/regionDir;
}


//- Read the states behind a validated header into a flat array
void readBody
(
    std::ifstream& is,
    scalarField& data,
    const fileName& statePath
)
{
    is.read
    (
        reinterpret_cast<char*>(data.data()),
        std::streamsize(data.size()*sizeof(scalar))
    );
    if (!is)
    {
        FatalErrorInFunction
            << "Unexpected end of restart state file " << statePath
            << exit(FatalError);
    }
    for (const scalar value : data)
    {
        restartStateIO::checkValue(value, statePath);
    }
}


//- Write a state file: header plus a flat array of rows*stride values
void writeState
(
    const fileName& statePath,
    const word& model,
    const label stride,
    const label rows,
    const scalarField& data
)
{
    mkDir(statePath.path());
    std::ofstream os(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }
    restartStateIO::writeHeader(model, stride, rows, os);
    os.write
    (
        reinterpret_cast<const char*>(data.cdata()),
        std::streamsize(data.size()*sizeof(scalar))
    );
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }
}


//- Open a state file and read its header
restartStateIO::Header openState
(
    const fileName& statePath,
    std::ifstream& is
)
{
    is.open(statePath.c_str(), std::ios::binary);
    if (!is)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath
            << exit(FatalError);
    }
    return restartStateIO::readHeader(is, statePath);
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Gather the per-processor restart state files (<model>State) of"
        " cardiacFoam's ionic and active-tension models into the undecomposed"
        " case, or scatter them into the processor directories with"
        " -decompose. decomposePar and reconstructPar leave these files alone."
    );
    timeSelector::addOptions(true, true);
    argList::addBoolOption
    (
        "decompose",
        "Scatter <time>/<model>State into the processor directories"
        " instead of gathering"
    );
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const word regionName =
        args.getOrDefault<word>("region", polyMesh::defaultRegion);
    const word regionDir = polyMesh::regionName(regionName);
    const bool decompose = args.found("decompose");

    const label nProcs =
    (
        regionDir.empty()
      ? fileHandler().nProcs(args.path())
      : fileHandler().nProcs(args.path(), regionName)
    );
    if (nProcs == 0)
    {
        FatalErrorInFunction
            << "No processor directories found in " << args.path()
            << exit(FatalError);
    }
    if (!regionDir.empty())
    {
        Info<< "Using region: " << regionName << nl;
    }
    Info<< "Processors: " << nProcs << nl << endl;

    // Processor databases and cell addressing
    PtrList<Time> procTimes(nProcs);
    List<labelList> addressing(nProcs);
    label nGlobalCells = 0;
    forAll(procTimes, proci)
    {
        procTimes.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/("processor" + Foam::name(proci))
            )
        );
        labelIOList addr
        (
            IOobject
            (
                "cellProcAddressing",
                procTimes[proci].constant(),
                regionDir/polyMesh::meshSubDir,
                procTimes[proci],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        addressing[proci].transfer(addr);
        if (addressing[proci].size())
        {
            nGlobalCells = max(nGlobalCells, max(addressing[proci]) + 1);
        }
    }

    label nProcCells = 0;
    forAll(addressing, proci)
    {
        nProcCells += addressing[proci].size();
    }
    if (nProcCells != nGlobalCells)
    {
        FatalErrorInFunction
            << "cellProcAddressing covers " << nProcCells
            << " cells but addresses " << nGlobalCells
            << exit(FatalError);
    }

    const instantList timeDirs =
    (
        decompose
      ? timeSelector::select(runTime.times(), args)
      : timeSelector::select(procTimes[0].times(), args)
    );
    if (timeDirs.empty())
    {
        Info<< "No times selected" << nl << endl;
        return 0;
    }

    label nFiles = 0;

    for (const instant& t : timeDirs)
    {
        const word& timeName = t.name();
        const fileName serialDir = regionTimeDir(args.path(), timeName, regionDir);

        if (decompose)
        {
            for (const word& name : stateFiles(serialDir))
            {
                const fileName statePath = serialDir/name;
                std::ifstream is;
                const restartStateIO::Header header = openState(statePath, is);
                if (header.rows != nGlobalCells)
                {
                    Info<< "    " << name << " at time " << timeName
                        << ": " << header.rows << " rows for " << nGlobalCells
                        << " cells, not a per-cell state, skipped" << nl;
                    continue;
                }
                scalarField data(header.rows*header.stride);
                readBody(is, data, statePath);

                forAll(addressing, proci)
                {
                    const labelList& addr = addressing[proci];
                    scalarField procData(addr.size()*header.stride);
                    forAll(addr, celli)
                    {
                        for (label s = 0; s < header.stride; ++s)
                        {
                            procData[celli*header.stride + s] =
                                data[addr[celli]*header.stride + s];
                        }
                    }
                    writeState
                    (
                        regionTimeDir(procTimes[proci].path(), timeName, regionDir)/name,
                        header.model,
                        header.stride,
                        addr.size(),
                        procData
                    );
                }
                Info<< "    " << name << " at time " << timeName << ": "
                    << header.rows << " cells x " << header.stride
                    << " states to " << nProcs << " processors" << nl;
                ++nFiles;
            }
        }
        else
        {
            const fileName proc0Dir =
                regionTimeDir(procTimes[0].path(), timeName, regionDir);

            for (const word& name : stateFiles(proc0Dir))
            {
                word model;
                label stride = -1;
                scalarField data;
                bool skip = false;

                forAll(addressing, proci)
                {
                    const labelList& addr = addressing[proci];
                    const fileName statePath =
                        regionTimeDir(procTimes[proci].path(), timeName, regionDir)/name;
                    std::ifstream is;
                    const restartStateIO::Header header = openState(statePath, is);

                    if (header.rows != addr.size())
                    {
                        Info<< "    " << name << " at time " << timeName
                            << ": processor" << proci << " has " << header.rows
                            << " rows for " << addr.size()
                            << " cells, not a per-cell state, skipped" << nl;
                        skip = true;
                        break;
                    }
                    if (proci == 0)
                    {
                        model = header.model;
                        stride = header.stride;
                        data.resize(nGlobalCells*stride);
                    }
                    else if (header.model != model || header.stride != stride)
                    {
                        FatalErrorInFunction
                            << statePath << " holds " << header.model << " x "
                            << header.stride << ", processor0 holds "
                            << model << " x " << stride
                            << exit(FatalError);
                    }

                    scalarField procData(addr.size()*stride);
                    readBody(is, procData, statePath);
                    forAll(addr, celli)
                    {
                        for (label s = 0; s < stride; ++s)
                        {
                            data[addr[celli]*stride + s] =
                                procData[celli*stride + s];
                        }
                    }
                }
                if (skip)
                {
                    continue;
                }

                writeState(serialDir/name, model, stride, nGlobalCells, data);
                Info<< "    " << name << " at time " << timeName << ": "
                    << nGlobalCells << " cells x " << stride << " states from "
                    << nProcs << " processors" << nl;
                ++nFiles;
            }
        }
    }

    if (nFiles == 0)
    {
        Info<< "No restart state files found in the selected times" << nl;
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
