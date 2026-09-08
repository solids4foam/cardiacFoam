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
    interpolateFibreField

Description
    Face-interpolates an existing cell-centred fibre direction field (f0)
    to a unit-length face field (f0f), required by electroMechanicalLaw.

    Unlike setFibreField, this does not compute f0 itself from a
    Laplace-solved transmural coordinate and a helix-angle rule: it reads
    whatever f0 is already present (e.g. supplied from an anatomical fibre
    dataset) and only produces its face-interpolated companion. Meant to
    be run once per mesh; commit the resulting f0f as source data rather
    than regenerating it every run.

Author
    Philip Cardiff, UCD (original setFibreField, whose interpolate+
    normalize pattern this mirrors).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading f0" << endl;
    volVectorField f0
    (
        IOobject
        (
            "f0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Interpolating f0 to faces and normalizing" << endl;
    surfaceVectorField f0f
    (
        IOobject
        (
            "f0f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(f0)
    );
    f0f /= mag(f0f);

    Info<< "Writing f0f" << endl;
    f0f.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
