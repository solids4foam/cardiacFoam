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
    Reads constant/<region>/polyMesh and detects whether the mesh points are in
    SI metres. Detect-only by default (no write). Pass -rescale to apply the
    auto-detected factor, or -scale <factor> to apply an explicit one; either
    rewrites the region to SI metres.

Usage
    checkMeshGeometry [-region <name>] [-rescale | -scale <factor>]

Author
    Simao Nieto de Castro. All rights reserved.
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
        "Detect mesh point units and rescale a region to SI metres."
    );
    argList::addBoolOption
    (
        "rescale",
        "Apply the auto-detected scale factor and rewrite the mesh "
        "(default: report only, no write)"
    );
    argList::addOption
    (
        "region",
        "name",
        "Mesh region to operate on (default: region0)"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Apply this explicit scale factor and rewrite (overrides auto-detect)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== checkMeshGeometry ==========\n" << endl;

    const word regionName =
        args.getOrDefault<word>("region", polyMesh::defaultRegion);

    polyMesh mesh
    (
        IOobject
        (
            regionName,
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

    // Unit-detection thresholds. MUST mirror
    // openfoam_driver/specs/mesh_geometry.py (guarded by
    // test_mesh_geometry_contract.py). mmLower is 20 (not 1) so large SI
    // domains (whole-torso meshes, ~1-2 m) are not mis-detected as mm.
    const scalar mmLower = 20.0;
    const scalar umLower = 1000.0;
    const scalar umUpper = 1e6;

    scalar scaleFactor = 1.0;
    word detectedUnit = "m";

    if (maxDim >= mmLower && maxDim < umLower)
    {
        scaleFactor = 1e-3;
        detectedUnit = "mm";
    }
    else if (maxDim >= umLower && maxDim < umUpper)
    {
        scaleFactor = 1e-6;
        detectedUnit = "um";
    }

    if (args.found("scale"))
    {
        scaleFactor = args.get<scalar>("scale");
        detectedUnit = "explicit";
    }

    // Detect-only is the default. A write requires explicit opt-in so the
    // utility never silently rescales a (possibly non-dimensional) mesh.
    const bool doWrite = args.found("scale") || args.found("rescale");

    if (mag(scaleFactor - 1.0) < SMALL)
    {
        Info<< "Mesh appears to be in meters (factor 1). Nothing to do." << nl
            << endl;
    }
    else if (!doWrite)
    {
        WarningInFunction
            << "Max dimension = " << maxDim
            << " suggests mesh is in " << detectedUnit << ", not meters.\n"
            << "  Detect-only (default): no write performed.\n"
            << "  Re-run with -rescale (auto) or -scale " << scaleFactor
            << " to apply." << nl << endl;
    }
    else
    {
        Info<< "Applying scale factor " << scaleFactor
            << " (" << detectedUnit << ") and rewriting constant/polyMesh ..."
            << nl << endl;

        pointField newPoints(mesh.points());
        newPoints *= scaleFactor;
        mesh.movePoints(newPoints);
        mesh.write();

        const boundBox bbScaled(mesh.points(), false);
        Info<< "Scaled bounding box : " << bbScaled << nl
            << "Scaled max dimension: " << bbScaled.maxDim() << " m"
            << nl << endl;
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
