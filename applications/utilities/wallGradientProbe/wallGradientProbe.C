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
    wallGradientProbe

Description
    Cell gradient of a prescribed field on a sealed-wall mesh through the
    sealedWallLeastSquares scheme, without a solve.

    The field is affine or the monodomain shell manufactured solution F.
    Sealed patches carry the adjacent cell value (the zero-gradient trace),
    every other non-coupled patch the exact value. Prints:
    - stockDiff: |sealedWallLeastSquares - leastSquares| with every patch
      stock (zero: the scheme is leastSquares there);
    - nameDiff: |grad(psi) - grad(psi.component(0))| on a field of another
      name (zero: both gradient names find the states of psi);
    - wall and interior RMS and maximum errors, relative to max |grad|.

    The case's gradSchemes select sealedWallLeastSquares for
    grad(psi), grad(psi.component(0)) and grad(stock), leastSquares for
    grad(reference).

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "bitSet.H"
#include "IOmanip.H"
#include "sealedWallPolicy.H"
#include "manufacturedShellMonodomainReference.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Sealed-wall cell gradient of an affine or manufactured field"
    );
    argList::addOption("field", "affine|shell", "prescribed field (shell)");
    argList::addOption
    (
        "sealed",
        "wordRes",
        "sealed patches (endocardium epicardium)"
    );
    argList::addOption("value", "wordRes", "value patches (none)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word field = args.getOrDefault<word>("field", "shell");
    if (field != "affine" && field != "shell")
    {
        FatalErrorInFunction
            << "-field takes affine or shell, not " << field
            << exit(FatalError);
    }
    const wordRes sealedNames
    (
        args.getOrDefault<wordRes>
        (
            "sealed",
            wordRes({wordRe("endocardium"), wordRe("epicardium")})
        )
    );
    const wordRes valueNames
    (
        args.getOrDefault<wordRes>("value", wordRes())
    );

    const shellMonodomainMMS::Params p = shellMonodomainMMS::defaultParams();
    const vector b(0.31, -0.47, 0.83);

    auto exact = [&](const point& x) -> scalar
    {
        return field == "affine" ? 0.2 + (b & x) : shellMonodomainMMS::F(x, p);
    };
    auto exactGrad = [&](const point& x) -> vector
    {
        return field == "affine" ? b : shellMonodomainMMS::gradF(x, p);
    };

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList states(patches.size(), label(sealedWallPolicy::stock));
    for (const label patchi : patches.indices(sealedNames))
    {
        states[patchi] = sealedWallPolicy::sealed;
    }
    for (const label patchi : patches.indices(valueNames))
    {
        states[patchi] = sealedWallPolicy::value;
    }
    sealedWallPolicy::publish(mesh, "psi", states);
    sealedWallPolicy::publish
    (
        mesh,
        "stock",
        labelList(patches.size(), label(sealedWallPolicy::stock))
    );

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::calculatedType()
    );
    forAll(psi, celli)
    {
        psi[celli] = exact(mesh.C()[celli]);
    }
    forAll(psi.boundaryField(), patchi)
    {
        fvPatchScalarField& psiP = psi.boundaryFieldRef()[patchi];
        if (psiP.coupled())
        {
            continue;
        }
        const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
        const labelUList& faceCells = psiP.patch().faceCells();
        forAll(psiP, i)
        {
            psiP[i] =
                states[patchi] == sealedWallPolicy::sealed
              ? psi[faceCells[i]]
              : exact(Cf[i]);
        }
    }
    psi.correctBoundaryConditions();

    volScalarField renamed("exactForGradient(psi)", psi);

    const volVectorField g(fvc::grad(psi, "grad(psi)"));
    const volVectorField gComponent
    (
        fvc::grad(renamed, "grad(psi.component(0))")
    );
    const volVectorField gStock(fvc::grad(psi, "grad(stock)"));
    const volVectorField gReference(fvc::grad(psi, "grad(reference)"));

    const scalar stockDiff =
        gMax(mag(gStock.primitiveField() - gReference.primitiveField())());
    const scalar nameDiff =
        gMax(mag(g.primitiveField() - gComponent.primitiveField())());

    bitSet isWall(mesh.nCells());
    forAll(states, patchi)
    {
        if (states[patchi] != sealedWallPolicy::stock)
        {
            isWall.set(patches[patchi].faceCells());
        }
    }

    scalar sumW = 0, volW = 0, maxW = 0;
    scalar sumI = 0, volI = 0, maxI = 0;
    scalar gradMax = 0;
    forAll(g, celli)
    {
        const vector ge(exactGrad(mesh.C()[celli]));
        const scalar e = mag(g[celli] - ge);
        const scalar V = mesh.V()[celli];
        gradMax = max(gradMax, mag(ge));
        if (isWall.test(celli))
        {
            sumW += V*sqr(e);
            volW += V;
            maxW = max(maxW, e);
        }
        else
        {
            sumI += V*sqr(e);
            volI += V;
            maxI = max(maxI, e);
        }
    }
    reduce(sumW, sumOp<scalar>());
    reduce(volW, sumOp<scalar>());
    reduce(maxW, maxOp<scalar>());
    reduce(sumI, sumOp<scalar>());
    reduce(volI, sumOp<scalar>());
    reduce(maxI, maxOp<scalar>());
    reduce(gradMax, maxOp<scalar>());

    const label nCells = returnReduce(mesh.nCells(), sumOp<label>());
    const scalar h = Foam::cbrt((volW + volI)/nCells);

    Info<< setprecision(10)
        << "field " << field
        << " cells " << nCells
        << " h " << h
        << " stockDiff " << stockDiff
        << " nameDiff " << nameDiff
        << " wallRMS " << Foam::sqrt(sumW/max(volW, VSMALL))/gradMax
        << " wallMax " << maxW/gradMax
        << " interiorRMS " << Foam::sqrt(sumI/max(volI, VSMALL))/gradMax
        << " interiorMax " << maxI/gradMax << endl;

    return 0;
}


// ************************************************************************* //
