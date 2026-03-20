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
    checkMeshGeometry

Description
    Reads constant/polyMesh and checks whether the mesh points are in SI
    meters. If not, warns the user and rescales the mesh to meters based
    on the detected order of magnitude (mm -> 1e-3, um -> 1e-6).

Usage
    checkMeshGeometry

Author
    cardiacFoam developers.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "boundBox.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Check mesh geometry units and auto-scale to SI meters if needed."
    );
    argList::addBoolOption
    (
        "noScale",
        "Suppress automatic scaling even if units appear incorrect"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== checkMeshGeometry ==========\n" << endl;

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const boundBox bb(mesh.points(), false);
    const scalar maxDim = bb.maxDim();

    Info<< "Bounding box : " << bb << nl
        << "Max dimension: " << maxDim << " [raw units]" << nl << endl;

    // Threshold-based unit detection
    scalar scaleFactor = 1.0;
    word detectedUnit = "m";

    if (maxDim >= 1.0 && maxDim < 1000.0)
    {
        scaleFactor = 1e-3;
        detectedUnit = "mm";
    }
    else if (maxDim >= 1000.0 && maxDim < 1e6)
    {
        scaleFactor = 1e-6;
        detectedUnit = "um";
    }

    const bool noScale = args.found("noScale");

    if (mag(scaleFactor - 1.0) < SMALL)
    {
        Info<< "Mesh appears to be in meters. No scaling applied." << nl
            << endl;
    }
    else
    {
        WarningInFunction
            << "Max dimension = " << maxDim
            << " suggests mesh is in " << detectedUnit
            << ", not meters.\n"
            << "  Applying scale factor: " << scaleFactor << "\n"
            << "  Rewriting constant/polyMesh ..." << nl
            << endl;

        if (noScale)
        {
            Info<< "  -noScale flag set: skipping write." << nl << endl;
        }
        else
        {
            pointField newPoints(mesh.points());
            newPoints *= scaleFactor;
            mesh.movePoints(newPoints);
            mesh.write();

            const boundBox bbScaled(mesh.points(), false);
            Info<< "Scaled bounding box : " << bbScaled << nl
                << "Scaled max dimension: " << bbScaled.maxDim() << " m"
                << nl << endl;
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
