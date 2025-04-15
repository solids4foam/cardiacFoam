/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    extractIdealisedVentricleResults

Description
    Extrac results for the idealised ventricle case, defined as problem 2 in
    Land et al., 2015, Verification of cardiac mechanics software: benchmark
    problems and solutions for testing active and passive material behaviour,
    Proc. R. Soc. A.47120150641, http://doi.org/10.1098/rspa.2015.0641.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "pointFields.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    argList::noParallel();

    // Get times list
    const instantList Times = runTime.times();

    // Go to the latest time step
    runTime.setTime(Times[Times.size() - 1], Times.size() - 1);
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Re-read the mesh, if needed
    // Assume total Lagrangian approach so the mesh does not move
    //mesh.readUpdate();

    // Read the displacement field
    const volVectorField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Read the point displacement field
    pointMesh pMesh(mesh);
    const pointVectorField pointD
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        pMesh
    );

    // Calculate the deformed mesh points
    const pointField deformedPoints(mesh.points() + pointD);

    // Find the apex position (most negative z) on the endocardial (inside) and
    // epicardial (outside) surfaces
    const label insidePatchID = mesh.boundaryMesh().findPatchID("inside");
    if (insidePatchID == -1)
    {
        FatalError
            << "Cannot find the 'inside' patch" << abort(FatalError);
    }
    const label outsidePatchID = mesh.boundaryMesh().findPatchID("outside");
    if (outsidePatchID == -1)
    {
        FatalError
            << "Cannot find the 'outside' patch" << abort(FatalError);
    }

    // Inside surface apex
    scalar insideApexZ = GREAT;
    const labelList& insideMeshPoints =
        mesh.boundaryMesh()[insidePatchID].meshPoints();
    forAll(insideMeshPoints, pI)
    {
        const label pointID = insideMeshPoints[pI];
        insideApexZ = min(insideApexZ, deformedPoints[pointID].z());
    }

    // Outside surface apex
    scalar outsideApexZ = GREAT;
    const labelList& outsideMeshPoints =
        mesh.boundaryMesh()[outsidePatchID].meshPoints();
    forAll(outsideMeshPoints, pI)
    {
        const label pointID = outsideMeshPoints[pI];
        outsideApexZ = min(outsideApexZ, deformedPoints[pointID].z());
    }

    Info<< nl
        << "Apex:" << nl
        << "    inside = " << insideApexZ << nl
        << "    outside = " << outsideApexZ << endl;

    // Define uniformly-sampled u field
    const label NUM_SAMPLES = 100;
    scalarField u(NUM_SAMPLES, 0);
    const scalar uMax = -Foam::acos(5.0/20.0);
    const scalar uMin = -constant::mathematical::pi;
    const scalar uStep = (uMax - uMin)/(u.size() - 1);
    for(int i = 0; i < u.size(); ++i)
    {
        u[i] = uMin + i*uStep;
    }
    // Info<< nl << "u = " << u << endl;

    // Define undeformed coordinates along the mid-line (t = 0.5)
    vectorField midLine(u.size(), vector::zero);
    const scalar rs = (0.007 + 0.010)/2.0;
    const scalar rl = (0.017 + 0.020)/2.0;
    forAll(midLine, pI)
    {
        midLine[pI] = vector(rs*Foam::sin(u[pI]), 0, rl*Foam::cos(u[pI]));
    }
    // Info<< nl << "midLine (undeformed) = " << midLine << endl;

    // Find the cells containing the midLine points
    meshSearch searchEngine(mesh);
    labelList midLineCellIDs(midLine.size(), -1);
    label seedCellID = -1;
    forAll(midLineCellIDs, I)
    {
        midLineCellIDs[I] = searchEngine.findCell(midLine[I], seedCellID);
        seedCellID = midLineCellIDs[I];
    }
    // Info<< "midLineCellIDs = " << midLineCellIDs << endl;
    if (min(midLineCellIDs) == -1)
    {
        FatalError
            << "Could not find cells on the mid line" << abort(FatalError);
    }

    // Determine the displacement of the midLine points
    // We have many options for this interpolation: let's extrapolate from the
    // cell centre using the displacement gradient and see if it's good enough.
    // Ideally, we would blend between values from neighbouring cells to avoid
    // step jumps, but let's keep it simple first.
    const volTensorField gradD(fvc::grad(D));
    const vectorField& DI = D;
    const vectorField& CI = mesh.C();
    const tensorField& gradDI = gradD;
    vectorField midLineDisp(midLineCellIDs.size(), vector::zero);
    vectorField midLineCoord(midLineCellIDs.size(), vector::zero);
    forAll(midLineCellIDs, cI)
    {
        const label cellID = midLineCellIDs[cI];
        const vector d = midLine[cI] - CI[cellID];
        midLineDisp[cI] = DI[cellID] + (d & gradDI[cellID]);
        midLineCoord[cI] = midLine[cI] + midLineDisp[cI];
    }
    // Info<< "midLineDisp = " << midLineDisp << endl;
    // Info<< "midLineCoord = " << midLineCoord << endl;

    // Write the midLine to a file
    OFstream outFile("midLineDeformed.txt");
    Info<< nl
        << "Writing the deformed midline coordinates to midLineDeformed.txt"
        << endl;
    outFile
        << "# x y z" << endl;
    forAll(midLineCoord, cI)
    {
        outFile
            << midLineCoord[cI].x() << " "
            << midLineCoord[cI].y() << " "
            << midLineCoord[cI].z() << endl;
    }

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
