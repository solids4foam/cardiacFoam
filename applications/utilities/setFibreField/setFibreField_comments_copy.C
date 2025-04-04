/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Application
    setFibreField

Description
    Set the initial fibre field for cardiac simulations.

    The utility currently sets the fibre field for the Land et al. (2015)
    ellipsoidal ventricle benchmarks.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "transform.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Read the transmural distance field
    // This will be calculated using the given boundary conditions
    // -A: Note that transmural distance is the shortest distance from the endo to epicardinal layers

    // -A: note that t here is openfoams version of cout/print/logging utility... the endl creates a new line at the end
    Info<< "Reading t" << endl;
    //-A: volScalarField is a class within openfoam that creates a scalar field for the centre of each volume block.
    // This part of the code is creating a scalar field for t.

    // -A: IO here is for creating an object in openfoam. t here is the name if the object.
    // runTime.timeName() is for specifying the current time domain. I.e "0" or "140" ect.
    // mesh associates the field with a computational mesh.
    // IOobject::MUST_READ means there must be a t file in the "0" directory to read from.
    // IOobject::AUTO_WRITE ensure that t is automatically saved when creating results.
    volScalarField t
    (
        IOobject
        (
            "t",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Solve diffusion equation to calculate the transmural distance (t)
    // Set fixedValue = 0.0 on the endocardium and fixedValue = 1.0 on the
    // epicardium, with zeroGradient on other patches
    // -A: logging that we are solving the diffusion equation for t. We are treating it as a heat propogation throughout the mesh to
    // get the values within the mesh.
    Info<< "Solving the diffusion equation for t" << endl;
    // -A: diffusion equation... steady state so just the laplacian? should this be steady state? because of incompressibility?
    // Does this make sense? for instance should it be steady state? or should the transmural directions evolve?
    fvScalarMatrix tEqn(fvm::laplacian(t));
    //-A: solving the diffusion equation.
    tEqn.solve();
    Info<< "Writing t" << endl;
    // -A: saving t values.
    t.write();

    // Base direction as per Rossi-Lassila et al. approach
    const vector k(0, 0, 1);

    // Transmural direction
    //-A: calculating the vector field of the gradient of transmural direction and normalising it
    //-A: does the below make sense as the transmural direction? Since t is defined initially as the transmural direction? or perhaps this is just used
    //-A: later to arrive at the transmural direction.
    volVectorField et("et", fvc::grad(t));
    et /= mag(et);

    // Normal direction -A: the scalar projection of transmural on k... i.e the part of transmural direction that lies on k.
    volVectorField en("en", k - (k & et)*et);
    en /= mag(en);

    // Longitudinal direction.... cross product.
    volVectorField el("el", en ^ et);
    el /= mag(el);

    // writing out the files to each one.
    Info<< "Writing et, en, el" << endl;
    et.write();
    en.write();
    el.write();

    //-A: perhaps the above is a decomposition of t into 3 components of the vector... i.e in each direction? but its missing one.

    // -A: this should be 90 to 90 I think?
    // alpha at the endocardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEndo = degToRad(60.0);

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEpi = degToRad(-60.0);

    // Calculate alphaRadians

    //-A: creating a scalar field called alphaRadians. Linearly interpolating between alpha epi and endo as f varies linearly from 0 to 1 here.
    // -A: I am not sure if this is supposed to be the alpha(t) from the paper which is defined to be 90 -180 t? or if this is strictly for the rotation
    // to get the mesh.
    const volScalarField alphaRadians
    (
        "alphaRadians",
        alphaEndo*(1.0 - t) + alphaEpi*t
    );
    Info<< "Writing alphaRadians" << endl;
    alphaRadians.write();

    // Rotate el about et by alphaRadians to get f
    // using alpha radians to get f? so I think it is the above?
    // I think alpha radians may be incorrect here since it should depend on t and is defined by 90 - 180t... but how can we use it to define t if
    // it is dependent on t? and i dont see a dependence on alpha endo or alpha api for alpha radians in the paper.
    // I am just trying to think this through here.. it could be correct but I need to make sure.
    // Perhaps this comes from the formula for f where we have f = n (dx/du)cos(alpha) + n (dx/dv) sin(alpha) since cos and sin vary from 1 to -1?


    //initialising the f0 to el field
    volVectorField f0("f0", el);
    //looping though the values of et for each cell.
    forAll(et, cellI)
    {
        const tensor rotT = Ra(et[cellI], alphaRadians[cellI]); //computing the rotation matrix for.. rotation around the longitudinal by alpharadians. 
        f0[cellI] = transform(-rotT, f0[cellI]); // Applying the rotation.
    }
    forAll(et.boundaryField(), patchI) // same idea as above but applied to the boundary.
    {
        forAll(et.boundaryField()[patchI], faceI)
        {
            const tensor rotT = Ra
            (
                et.boundaryField()[patchI][faceI],
                alphaRadians.boundaryField()[patchI][faceI]
            );
            f0.boundaryFieldRef()[patchI][faceI] = transform
            (
                -rotT,
                f0.boundaryField()[patchI][faceI]
            );
        }
    }

    f0.write();


    // Calculate surface fields
    const surfaceScalarField tf(fvc::interpolate(t));

    // Transmural direction
    surfaceVectorField etf("etf", fvc::interpolate(fvc::grad(t)));
    //surfaceVectorField n(mesh.Sf()/mesh.magSf());
    // etf += fvc::snGrad(t)*n - ((I - sqr(n)) & etf);
    etf /= mag(etf);

    // Normal direction
    surfaceVectorField enf("enf", k - (k & etf)*etf);
    enf /= mag(enf);

    // Longitudinal direction
    surfaceVectorField elf("elf", enf ^ etf);
    elf /= mag(elf);

    etf.write();
    enf.write();
    elf.write();

    // Calculate alphaRadians
    const surfaceScalarField alphaRadiansf
    (
        "alphaRadiansf",
        alphaEndo*(1.0 - tf) + alphaEpi*tf
    );
    Info<< "Writing alphaRadiansf" << endl;
    alphaRadiansf.write();

    // Rotate el about et by alphaRadians to get f
    surfaceVectorField f0f("f0f", elf);
    forAll(etf, faceI)
    {
        const tensor rotT = Ra(etf[faceI], alphaRadiansf[faceI]);
        f0f[faceI] = transform(-rotT, f0f[faceI]);
    }
    forAll(etf.boundaryField(), patchI)
    {
        forAll(etf.boundaryField()[patchI], faceI)
        {
            const tensor rotT = Ra
            (
                etf.boundaryField()[patchI][faceI],
                alphaRadiansf.boundaryField()[patchI][faceI]
            );
            f0f.boundaryFieldRef()[patchI][faceI] = transform
            (
                -rotT,
                f0f.boundaryField()[patchI][faceI]
            );
        }
    }

    f0f.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
