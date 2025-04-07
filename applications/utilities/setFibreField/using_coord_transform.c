/*---------A---A------A---------------------------------------------------------*\
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
    et /= (magSqr(et)); // transforming this into dx/dt.
    et /= mag(et);

    // Normal direction -A: the scalar projection of transmural on k... i.e the part of transmural direction that lies on k.
    volVectorField en("en", k - (k & et)*et);
    en /= mag(en);

    // Longitudinal direction.... cross product.
    volVectorField el("el", en ^ et);
    el /= mag(el);

    // writing out the files to each one.
    Info<< "Writing et, en, el, ec" << endl;
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
	//degToRadian(90 - 180*t)
	//(degToRad(90) - degToRad(180) * t)
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


    /*
     *
     *TRYING A NEW APPROACH BY EXPLICITLY FOLLOWING THE NEWEST PAPER.
     *
     *
     */

    //creating rs and rs to define other stuff.

    // Eq (16) from the paper
    // with r_long_endo = 9x10^-2
    // with r_long_epi = 9.7x10^-2
    // with r_short_endo = 2.5x10^-2
    // with r_short_epi = 3.5x10^-2 
    const volScalarField rs
    (
        "rs",
	(2.5 + 1*t) // * 0.01
    );
    Info<< "Writing rs" << endl;
    rs.write();

    
    const volScalarField rl
    (
        "rl",
        (9 + 0.7* t) // * 0.01
    );
    Info<< "Writing rl" << endl;
    rl.write();

    
    //initialising mu from the paper
    volScalarField mu
      (
       IOobject
       (
        "mu",                    // Name of the field
        runTime.timeName(),     // Time directory (e.g., "0", "1", etc.)
        mesh,                   // Mesh to associate the field with
        IOobject::NO_READ,      // Do not read from disk (no existing file)
        IOobject::AUTO_WRITE    // Write the field automatically to disk
	),
       mesh,
       dimensionedScalar("mu", dimless, 0.0) // Mesh to define the field on
       );

    //Initialising thetha from the paper
    volScalarField thetha
      (
       IOobject
       (
        "thetha",                    // Name of the field                                                                                                        
        runTime.timeName(),     // Time directory (e.g., "0", "1", etc.)                                                                                    
        mesh,                   // Mesh to associate the field with                                                                                         
        IOobject::NO_READ,      // Do not read from disk (no existing file)                                                                                 
        IOobject::AUTO_WRITE    // Write the field automatically to disk                                                                                    
        ),
       mesh,
       dimensionedScalar("thetha", dimless, 0.0)                                                                                             
       );

    //looping through an updating the mu and thetha values based procedure described in paper
    forAll(mesh.C(), cellI)
      {
	// Here I was thinking that x, y,z are supposed to come from the local
	// coordinate system.
	//vector m = ((mesh.C()[cellI] & en[cellI]) * en[cellI])
	  // +((mesh.C()[cellI] & el[cellI]) * el[cellI]);
	
	scalar x =  el[cellI].x() + en[cellI].x();
	scalar y = el[cellI].y() + en[cellI].y();
	scalar z = el[cellI].z() + el[cellI].z();
	
	scalar muValue, thethaValue;
	scalar a, b;
	a = Foam::sqrt((y * y) + (z * z)) / (rs[cellI]); //parameter given in paper to 
	b = x / (rl[cellI]);                             // get mu and thetha 
	muValue = Foam::atan2(a, b);
	mu[cellI] = -muValue;

	if(muValue <= 1e-7)
	  {
	    thethaValue = 0.0;
	    thetha[cellI] = thethaValue;
	  }
	
	else
	  {
	    thethaValue = Foam::constant::mathematical::pi - Foam::atan2(z, - y);
	    thetha[cellI] = thethaValue;
	  }
      }

    forAll(mu.boundaryField(), patchI)  // Loop over boundary patches
      {
	forAll(mu.boundaryField()[patchI], faceI)  // Loop over faces in patch
	  {
	    //I tried to do the same x,y,z from the local coord but got strange results
	    scalar x =  el.boundaryField()[patchI][faceI].x() + en.boundaryField()[patchI][faceI].x();
	    scalar y = el.boundaryField()[patchI][faceI].y() + el.boundaryField()[patchI][faceI].y();
	    scalar z = el.boundaryField()[patchI][faceI].z() + el.boundaryField()[patchI][faceI].z();
	    scalar a = Foam::sqrt((y * y) + (z * z)) / (rs.boundaryField()[patchI][faceI]);
	    scalar b = x / (rl.boundaryField()[patchI][faceI]);

	    
	    scalar muValue = Foam::atan2(a, b);
	    mu.boundaryFieldRef()[patchI][faceI] = -muValue;

	    scalar thethaValue;

	    if (muValue <= 1e-7)
	      {
		thethaValue = 0.0;
		thetha.boundaryFieldRef()[patchI][faceI] = thethaValue;
	      }
	    else
	      {
		thethaValue = Foam::constant::mathematical::pi -  Foam::atan2(z, - y);
		thetha.boundaryFieldRef()[patchI][faceI] = thethaValue;
	      }
	  }
      }

    
    Info<< "Writing mu" << endl;
    mu.write();
    Info<< "Writing thetha" << endl;
    thetha.write();

    
    volVectorField e_mu("e_mu", el);

    forAll(mesh.C(), cellI) {
      //NOTE I have changed the x and z values for the derivative since it
      //gave more reasonable results.

      //Below is the derivative of each component of the vector from
      //Eq (12) w.r.t mu.
      scalar e_mu_x = rs[cellI] * Foam::cos(mu[cellI]) * Foam::sin(thetha[cellI]);
      scalar e_mu_y = rs[cellI] * Foam::cos(mu[cellI]) * Foam::cos(thetha[cellI]);
      scalar e_mu_z = -rl[cellI] * Foam::sin(mu[cellI]);

      // Assign the value to the dx_du field (as a vector)
      e_mu[cellI] = vector(e_mu_x, e_mu_y, e_mu_z); // Create a vector from the components
    }

    //Doing the same procedure for the boundaries
    forAll(e_mu.boundaryField(), patchI)
      {
	forAll(e_mu.boundaryField()[patchI], faceI)
	  {
	    scalar e_mu_x = rs.boundaryField()[patchI][faceI] *
	                     Foam::cos(mu.boundaryField()[patchI][faceI]) *
	                     Foam::sin(thetha.boundaryField()[patchI][faceI]);

	    scalar e_mu_y = rs.boundaryField()[patchI][faceI]
	                     *Foam::cos(mu.boundaryField()[patchI][faceI])
	                     *Foam::cos(thetha.boundaryField()[patchI][faceI]);

	    scalar e_mu_z = -rl.boundaryField()[patchI][faceI]
	                     *Foam::sin(mu.boundaryField()[patchI][faceI]);

	    // Assign the computed vector to the boundary field
	    e_mu.boundaryFieldRef()[patchI][faceI] = vector(e_mu_x, e_mu_y, e_mu_z);
	  }
      }

    e_mu /= mag(e_mu);
    Info<< "Writing e_mu" << endl;
    e_mu.write();

    

    volVectorField e_thetha("e_thetha", el);

    //NOTE I HAVE CHANGED X AND Z FOR THE FOLLOWING AS WELL.
    forAll(mesh.C() ,cellI) {

      // The following is the derivative for the x,y,z component w.r.t thetha in Eq (12) 
      scalar e_thetha_x = rs[cellI] * Foam::sin(mu[cellI]) * Foam::cos(thetha[cellI]);
      scalar e_thetha_y = -rs[cellI] * Foam::sin(mu[cellI]) * Foam::sin(thetha[cellI]);
      scalar e_thetha_z = 0;           // 0

      e_thetha[cellI] = vector(e_thetha_x, e_thetha_y, e_thetha_z); // Create a vector from the components
    }

    //Doing the same procedure for the boundaries
    forAll(e_thetha.boundaryField(), patchI)
      {
	forAll(e_thetha.boundaryField()[patchI], faceI)
	  {
	    scalar e_thetha_x = rs.boundaryField()[patchI][faceI]
	                      * Foam::sin(mu.boundaryField()[patchI][faceI])
	                      * Foam::cos(thetha.boundaryField()[patchI][faceI]);

	    scalar e_thetha_y = -rs.boundaryField()[patchI][faceI]
	                       *Foam::sin(mu.boundaryField()[patchI][faceI])
	                       *Foam::sin(thetha.boundaryField()[patchI][faceI]);

	    scalar e_thetha_z = 0.0; // No change in the z-direction

	    // Assign the computed vector to the boundary field
	    e_thetha.boundaryFieldRef()[patchI][faceI] = vector(e_thetha_x, e_thetha_y, e_thetha_z);
	  }
      }

    e_thetha /= mag(e_thetha);
    Info<< "Writing e_thetha" << endl;
    e_thetha.write();

    //Initialising the fibres
    volVectorField f1("f1", el);

    //Calculating the fibres based on Eq (15)
    // also projecting these results onto local coordinate system
    // well attempting to at least.
    forAll(mesh.C(), cellI)
    {
      f1[cellI] = (
		   ((( e_mu[cellI] * Foam::sin(alphaRadians[cellI]))
		     + ( e_thetha[cellI] * Foam::cos(alphaRadians[cellI]))) & el[cellI] )*el[cellI]
		   )
	
	         + (
		    ((( e_mu[cellI] * Foam::sin(alphaRadians[cellI])) +
		      ( e_thetha[cellI] * Foam::cos(alphaRadians[cellI]))) & en[cellI]) * en[cellI]
		    );

    }

    // doing the same procedure for the boundaries
    forAll(f1.boundaryField(), patchI)
      {
	forAll(f1.boundaryField()[patchI], faceI)
	  {
	    f1.boundaryFieldRef()[patchI][faceI] = (
						    (
						     ((e_mu.boundaryField()[patchI][faceI]
						       *Foam::sin(alphaRadians.boundaryField()[patchI][faceI])) + 
						      (e_thetha.boundaryField()[patchI][faceI] *
						       Foam::cos(alphaRadians.boundaryField()[patchI][faceI])))&
						     el.boundaryField()[patchI][faceI]
						     ) * el.boundaryField()[patchI][faceI]
						    )
	      
	                                         +
	      (
	       (
		((e_mu.boundaryField()[patchI][faceI]
		  *Foam::sin(alphaRadians.boundaryField()[patchI][faceI])) +
		  (e_thetha.boundaryField()[patchI][faceI] *
		   Foam::cos(alphaRadians.boundaryField()[patchI][faceI])))
		 & en.boundaryField()[patchI][faceI]
		 ) * en.boundaryField()[patchI][faceI]
	       );
	  }
	  }

    // Optionally, write out the result to a file
    f1 /= mag(f1);
    Info<< "Writing f1" << endl;
    f1.write();


    //writing out the sheet normal and the sheet directions.
    //From Eq (15) with an extra calculation for the projection of coordinated.
    volVectorField n_("n_", (( ( e_mu & en )* en) + ( (e_mu & el ) * el ) ) ^ ( (( e_thetha & en )* en) + ( (e_thetha & el ) * el )) );
    n_ /= mag(n_);

    volVectorField s_("s_", f1 ^ n_ );
    s_ /= mag(s_);

    Info<< "Writting n_ and s_" << endl;
    n_.write();
    s_.write();





    

    /*End of trying new approach */












    
    
    //initialising the f0 to el field
    volVectorField f0("f0", f1);
    //looping though the values of et for each cell.
    /*
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
    */
    f0 /= mag(f0);
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
