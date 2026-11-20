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
    pvjContinuumParity

Description
    Evaluates pvjContinuumReference (src/verificationModels/coupledVerification)
    on the mesh of a case, so that the NumPy module
    mms_common/pvj_continuum.py (emitted as reference/pvj_continuum.py) can be
    compared with the C++ it ports. The class is compiled into the application;
    no solver library is needed.

Usage
    pvjContinuumParity -case <case>

    The case holds a polyMesh and system/pvjParityDict:

        kernel      uniform | linear | gaussian | cosine;
        radius      R;
        terminal    (x y z);
        field       none | cosine | poly;   // f of the cell and support integrals
        support
        {
            kind        interior | plane | cavitySphere;
            normal      (nx ny nz);          // plane, cavitySphere
            wallRadius  rho;                 // cavitySphere
        }

        // optional, as read by pvjContinuumReference
        referenceRadialPoints, referencePolarPoints, referenceAzimuthalPoints,
        referenceCutSize, referenceInsideSize, referenceTipSize

    Output, one record per line, 17 significant digits; other lines (the
    OpenFOAM banner) are to be ignored:

        G <cell> <volume> <Cx> <Cy> <Cz>     every cell: mesh.V(), mesh.C()
        S <W> <int w f over the support>     supportIntegral (f = 1 / field)
        I <cell> <omega> <omegaF>            cellIntegrals, cells that meet B

    field cosine = cos(pi x) cos(2 pi y) cos(3 pi z) (the tissue field of the
    cube junction at t = 0); poly = 1 + 0.3x - 0.7y + 0.2z + xy - 0.5z^2.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"

#include "pvjContinuumReference.H"

#include <cstdio>

using namespace Foam;

int main(int argc, char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    IOdictionary dict
    (
        IOobject
        (
            "pvjParityDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word kernel(dict.get<word>("kernel"));
    const scalar R = dict.get<scalar>("radius");
    const point terminal(dict.get<point>("terminal"));
    const word field(dict.getOrDefault<word>("field", "none"));
    const dictionary& supportDict = dict.subDict("support");

    const scalar pi = constant::mathematical::pi;
    pvjContinuumReference::scalarFunction F;
    if (field == "cosine")
    {
        F = [pi](const point& x)
        {
            return Foam::cos(pi*x.x())*Foam::cos(2*pi*x.y())*Foam::cos(3*pi*x.z());
        };
    }
    else if (field == "poly")
    {
        F = [](const point& x)
        {
            return 1 + 0.3*x.x() - 0.7*x.y() + 0.2*x.z() + x.x()*x.y()
                 - 0.5*x.z()*x.z();
        };
    }
    else if (field != "none")
    {
        FatalErrorInFunction
            << "field must be none, cosine or poly; found " << field
            << exit(FatalError);
    }

    const pvjContinuumReference ref(kernel, R, dict);

    const vectorField& centres = mesh.C();
    const scalarField& volumes = mesh.V();
    forAll(centres, cellI)
    {
        const point& C = centres[cellI];
        std::printf
        (
            "G %d %.17g %.17g %.17g %.17g\n",
            int(cellI), volumes[cellI], C.x(), C.y(), C.z()
        );
    }

    pvjContinuumReference::support s;
    s.centre = terminal;
    s.normal = supportDict.getOrDefault<vector>("normal", vector::zero);
    s.wallRadius = supportDict.getOrDefault<scalar>("wallRadius", 0);
    const word kind(supportDict.get<word>("kind"));
    if (kind == "interior")
    {
        s.kind = pvjContinuumReference::supportKind::interior;
    }
    else if (kind == "plane")
    {
        s.kind = pvjContinuumReference::supportKind::plane;
    }
    else if (kind == "cavitySphere")
    {
        s.kind = pvjContinuumReference::supportKind::cavitySphere;
    }
    else
    {
        FatalErrorInFunction
            << "support kind must be interior, plane or cavitySphere; found "
            << kind << exit(FatalError);
    }

    std::printf
    (
        "S %.17g %.17g\n",
        ref.supportIntegral(s),
        F ? ref.supportIntegral(s, F) : scalar(0)
    );

    labelList cells;
    scalarField omega, omegaF;
    ref.cellIntegrals(mesh, terminal, F, cells, omega, omegaF);

    forAll(cells, i)
    {
        std::printf
        (
            "I %d %.17g %.17g\n", int(cells[i]), omega[i], omegaF[i]
        );
    }

    return 0;
}


// ************************************************************************* //
