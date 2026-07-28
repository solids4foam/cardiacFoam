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
    const scalarField e(mag(gradNumeric.primitiveField() - gradExact.primitiveField()));

    scalar sumSq = 0.0;
    forAll(e, cellI)
    {
        sumSq += e[cellI]*e[cellI];
    }
    label nCells = e.size();

    reduce(sumSq, sumOp<scalar>());
    reduce(nCells, sumOp<label>());

    const scalar rms = nCells > 0 ? Foam::sqrt(sumSq/nCells) : 0.0;
    const scalar linf = gMax(e);

    // Reference magnitude, so the reported error can be read as relative.
    const scalar refRms =
        Foam::sqrt(gSum(magSqr(gradExact.primitiveField()))/max(nCells, 1));

    Info<< "gradient error against k*exp(k.x)" << nl
        << "  E_RMS,cell : " << rms << nl
        << "  E_inf      : " << linf << nl
        << "  |grad|_RMS : " << refRms << nl
        << "  relative   : " << (refRms > SMALL ? rms/refRms : 0.0) << nl
        << endl;

    Info<< "CSV," << nCells << "," << rms << "," << linf << "," << refRms << endl;

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
