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
    gradientReconstructionOrder

Description
    Measures the convergence of the cell-gradient reconstruction alone, with
    no solve involved.

    The eikonal manufactured solution is tau = exp(k.x), whose gradient is
    k*exp(k.x) in closed form. This utility sets tau exactly in the cells and
    on the boundary faces, applies fvc::grad through the scheme named in
    system/fvSchemes, and reports the error of the reconstructed gradient
    against the analytic one.

    The point is to separate two things the solved eikonal result cannot. A
    solved activation-time field carries the gradient error, the nonlinear
    fixed-point error and the boundary treatment together. Here the input
    field is exact by construction, so whatever error appears is the
    reconstruction operator and nothing else. Running this over the same mesh
    ladder as the eikonal study therefore tests directly whether the gradient
    operator is the ceiling on the activation-time order.

    Usage:
        gradientReconstructionOrder [-gradName <name>]

    -gradName selects the entry in fvSchemes/gradSchemes; default is
    "grad(tauExact)", falling back to the gradSchemes default.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "manufacturedEikonalReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addOption
    (
        "gradName",
        "word",
        "gradSchemes entry to use (default grad(tauExact))"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word gradName =
        args.getOrDefault<word>("gradName", "grad(tauExact)");

    // The manufactured wave vector, taken from the same reference header the
    // verifier uses so the two cannot drift apart.
    const vector k = manufacturedEikonalK(3);

    Info<< "Manufactured eikonal gradient reconstruction test" << nl
        << "  cells      : " << returnReduce(mesh.nCells(), sumOp<label>()) << nl
        << "  k          : " << k << nl
        << "  gradScheme : " << gradName << nl << endl;

    volScalarField tauExact
    (
        IOobject("tauExact", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless, Zero),
        "calculated"
    );

    // Exact cell values.
    const volVectorField& C = mesh.C();
    forAll(tauExact, cellI)
    {
        tauExact[cellI] = manufacturedEikonalTau(C[cellI], k);
    }

    // Exact boundary-face values, so that the only approximation left in the
    // reconstruction is the interior stencil rather than the boundary data.
    volScalarField::Boundary& tauBf = tauExact.boundaryFieldRef();
    forAll(tauBf, patchI)
    {
        fvPatchScalarField& pf = tauBf[patchI];

        if (pf.empty() || pf.type() == "empty")
        {
            continue;
        }

        const vectorField& Cf = pf.patch().Cf();
        forAll(pf, faceI)
        {
            pf[faceI] = manufacturedEikonalTau(Cf[faceI], k);
        }
    }

    // Reconstructed gradient through the production operator.
    const volVectorField gradNumeric
    (
        "gradNumeric",
        fv::gradScheme<scalar>::New
        (
            mesh,
            mesh.gradScheme(gradName)
        )().grad(tauExact, gradName)
    );

    // Analytic gradient of exp(k.x) is k*exp(k.x).
    volVectorField gradExact
    (
        IOobject("gradExact", runTime.timeName(), mesh),
        mesh,
        dimensionedVector(dimless, Zero),
        "calculated"
    );
    forAll(gradExact, cellI)
    {
        gradExact[cellI] = k*manufacturedEikonalTau(C[cellI], k);
    }

    // Cell-RMS and maximum norms of the vector error magnitude, matching the
    // norm definitions the verification models use.

    // We need to identify boundary cells
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    boolList isBoundaryCell(mesh.nCells(), false);
    forAll(bMesh, patchI)
    {
        const polyPatch& pp = bMesh[patchI];
        if (!pp.empty() && pp.type() != "empty")
        {
            const labelUList& faceCells = pp.faceCells();
            forAll(faceCells, i)
            {
                isBoundaryCell[faceCells[i]] = true;
            }
        }
    }

    const scalarField errMag(mag(gradNumeric.primitiveField() - gradExact.primitiveField()));
    const volScalarField::Internal& V = mesh.V();

    scalar errL2Vol = 0.0;
    scalar refL2Vol = 0.0;
    scalar errL2VolBulk = 0.0;
    scalar refL2VolBulk = 0.0;
    scalar errL2VolBound = 0.0;
    scalar refL2VolBound = 0.0;
    scalar errInf = 0.0;
    scalar errL1Vol = 0.0;
    scalar volTotal = 0.0;
    label cellsAbove = 0;

    forAll(errMag, cellI)
    {
        const scalar ev = errMag[cellI];
        const scalar rv = mag(gradExact[cellI]);
        const scalar v = V[cellI];

        if (ev > errInf) errInf = ev;
        if (ev > 0.05) cellsAbove++;

        errL1Vol += ev * v;
        const scalar ev2v = ev * ev * v;
        const scalar rv2v = rv * rv * v;

        errL2Vol += ev2v;
        refL2Vol += rv2v;
        volTotal += v;

        if (isBoundaryCell[cellI])
        {
            errL2VolBound += ev2v;
            refL2VolBound += rv2v;
        }
        else
        {
            errL2VolBulk += ev2v;
            refL2VolBulk += rv2v;
        }
    }

    reduce(errL2Vol, sumOp<scalar>());
    reduce(refL2Vol, sumOp<scalar>());
    reduce(errL2VolBulk, sumOp<scalar>());
    reduce(refL2VolBulk, sumOp<scalar>());
    reduce(errL2VolBound, sumOp<scalar>());
    reduce(refL2VolBound, sumOp<scalar>());
    reduce(errInf, maxOp<scalar>());
    reduce(errL1Vol, sumOp<scalar>());
    reduce(volTotal, sumOp<scalar>());
    reduce(cellsAbove, sumOp<label>());

    label nCells = errMag.size();
    reduce(nCells, sumOp<label>());

    const scalar rmsVol = volTotal > SMALL ? Foam::sqrt(errL2Vol/volTotal) : 0.0;
    const scalar rmsVolBulk = volTotal > SMALL ? Foam::sqrt(errL2VolBulk/volTotal) : 0.0;
    const scalar rmsVolBound = volTotal > SMALL ? Foam::sqrt(errL2VolBound/volTotal) : 0.0;
    const scalar meanE = volTotal > SMALL ? errL1Vol/volTotal : 0.0;
    const scalar refRmsVol = volTotal > SMALL ? Foam::sqrt(refL2Vol/volTotal) : 0.0;

    Info<< "gradient error against k*exp(k.x)" << nl
        << "  nCells = " << nCells << nl
        << "  E_inf = " << errInf << nl
        << "  Mean E = " << meanE << nl
        << "  Cells with error > 0.05 = " << cellsAbove << nl
        << "  L2 Bulk = " << rmsVolBulk << nl
        << "  L2 Bound = " << rmsVolBound << nl
        << "  L2 Total = " << rmsVol << nl
        << "  |grad|_RMS = " << refRmsVol << nl
        << "  relative = " << (refRmsVol > SMALL ? rmsVol/refRmsVol : 0.0) << nl
        << endl;

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
