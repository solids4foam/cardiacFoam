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
