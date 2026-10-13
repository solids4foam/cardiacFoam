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

#include "ecgModelIO.H"
#include "OSspecific.H"
#include <fstream>
#include <sstream>
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

autoPtr<OFstream> ecgModelIO::openTimeSeries
(
    const fileName& outDir,
    const word& filename,
    const wordList& electrodeNames
)
{
    mkDir(outDir);
    autoPtr<OFstream> osPtr(new OFstream(outDir/filename));

    OFstream& os = osPtr.ref();
    os.setf(std::ios::scientific);
    os.precision(8);

    os << "# time";
    forAll(electrodeNames, i)
    {
        os << "  " << electrodeNames[i];
    }
    os << nl;

    return osPtr;
}


autoPtr<OFstream> ecgModelIO::openRestartTimeSeries
(
    const fileName& outDir, const word& filename,
    const wordList& channelNames, scalar restartTime
)
{
    const fileName path(outDir/filename);
    if (restartTime <= 0 || !isFile(path))
        return openTimeSeries(outDir,filename,channelNames);
    std::ifstream previous(path.c_str());
    std::string line, token;
    std::getline(previous,line);
    std::istringstream header(line);
    header >> token;
    bool valid=token=="#";
    header >> token;
    valid=valid && token=="time";
    for (const word& name:channelNames)
    {
        if (!(header >> token) || token != name.c_str()) valid=false;
    }
    if (header >> token) valid=false;
    scalar lastTime=-GREAT;
    while (std::getline(previous,line))
    {
        if (line.empty()) continue;
        std::istringstream row(line);
        scalar time, value;
        if (!(row >> time) || !std::isfinite(time) || time<=lastTime)
        { valid=false; break; }
        forAll(channelNames,i) if (!(row >> value) || !std::isfinite(value)) valid=false;
        if (row >> token) valid=false;
        lastTime=time;
    }
    if (!valid || lastTime > restartTime+SMALL)
    {
        FatalErrorInFunction << "Cannot append " << filename
            << ": incompatible trace or existing samples later than the restart time."
            << exit(FatalError);
    }
    autoPtr<OFstream> out(new OFstream(path,IOstreamOption(),IOstreamOption::APPEND));
    out->setf(std::ios::scientific); out->precision(8);
    return out;
}


void ecgModelIO::writeRow
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


autoPtr<triSurface> ecgModelIO::loadSurface
(
    const dictionary& dict,
    const Time& runTime,
    pointField& faceCentres,
    fileName& vtkDir
)
{
    if (!dict.found("torsoSurface"))
    {
        return autoPtr<triSurface>(nullptr);
    }

    const fileName stlPath(dict.get<fileName>("torsoSurface"));
    const fileName fullPath(runTime.path()/stlPath);

    Info<< "ecgModelIO: loading torso surface from '"
        << fullPath << "'" << nl << endl;

    autoPtr<triSurface> surfPtr(new triSurface(fullPath));
    faceCentres = surfPtr().faceCentres();

    Info<< "Torso surface: " << faceCentres.size()
        << " faces loaded." << nl << endl;

    vtkDir = runTime.path()/"postProcessing"/"ECG";
    mkDir(vtkDir);

    return surfPtr;
}


void ecgModelIO::writeSurfaceVtk
(
    const triSurface& surf,
    const List<scalar>& phi,
    const fileName& vtkDir,
    const word& timeName,
    const word& fieldName
)
{
    const pointField& pts = surf.points();
    const label nPoints   = pts.size();
    const label nTris     = surf.size();

    const fileName vtkFile(vtkDir/fieldName + "_" + timeName + ".vtk");

    OFstream os(vtkFile);
    os.setf(std::ios::scientific);
    os.precision(8);

    os << "# vtk DataFile Version 3.0\n";
    os << "ECG potential (" << fieldName << ") t=" << timeName << "\n";
    os << "ASCII\n";
    os << "DATASET POLYDATA\n";

    os << "POINTS " << nPoints << " float\n";
    forAll(pts, pI)
    {
        os << pts[pI].x() << " " << pts[pI].y() << " " << pts[pI].z() << "\n";
    }

    os << "POLYGONS " << nTris << " " << 4*nTris << "\n";
    forAll(surf, tI)
    {
        const triFace& tri = surf[tI];
        os << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    os << "CELL_DATA " << nTris << "\n";
    os << "SCALARS " << fieldName << " float 1\n";
    os << "LOOKUP_TABLE default\n";
    forAll(phi, fI)
    {
        os << phi[fI] << "\n";
    }

    Info<< "ecgModelIO: surface written to " << vtkFile << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
