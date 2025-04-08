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
    Info<< "Reading t" << endl;
    
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

    Info<< "Solving the diffusion equation for t" << endl;
    
    fvScalarMatrix tEqn(fvm::laplacian(t));
 
    tEqn.solve();
    Info<< "Writing t" << endl;
  
    t.write();

    //using alpha radians as described in Arostica
    const scalar alphaEndo = 60.0;

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEpi = -60.0;

    const scalar pi = Foam::constant::mathematical::pi;
    // Calculate alphaRadians Eq (16) of Arostica
    const volScalarField alphaRadians
    (
        "alphaRadians",
        (alphaEndo + (alphaEpi - alphaEndo)*t) * pi/180
    );
    Info<< "Writing alphaRadians" << endl;
    alphaRadians.write();


    /*
     *
     *TRYING A NEW APPROACH BY EXPLICITLY FOLLOWING THE NEWEST PAPER.
     *
     *
     */

    //creating rs and rs to define other stuff.

    //Constants as defined from Eq's (10-11)
    // with r_long_endo = 9x10^-2
    const scalar r_long_endo = 9*0.01;
    // with r_long_epi = 9.7x10^-2
    const scalar r_long_epi = 9.7*0.01;
    // with r_short_endo = 2.5x10^-2
    const scalar r_short_endo = 2.5*0.01;
    // with r_short_epi = 3.5x10^-2
    const scalar r_short_epi = 3.5*0.01;

    //from Eq (16) of the paper we get rs and rl
    const volScalarField rs
    (
        "rs",
	r_short_endo + (r_short_epi - r_short_endo)*t
    );
    Info<< "Writing rs" << endl;
    rs.write();

    
    const volScalarField rl
    (
        "rl",
        r_long_endo + (r_long_epi - r_long_endo)*t
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
        
	scalar x =  mesh.C()[cellI].x();
	scalar y = mesh.C()[cellI].y();
	scalar z = mesh.C()[cellI].z();
	
	scalar muValue, thethaValue;
	scalar a, b;
	a = Foam::sqrt((y * y) + (z * z)) / (rs[cellI]); //parameter given in paper to 
	b = x / (rl[cellI]);                             // get mu and thetha 
	muValue = Foam::atan2(a, b);
	mu[cellI] = muValue;

	if(muValue <= 1e-7)
	  {
	    thethaValue = 0.0;
	    thetha[cellI] = thethaValue;
	  }
	
	else
	  {
	    thethaValue = pi - Foam::atan2(z, - y);
	    thetha[cellI] = thethaValue;
	  }
      }

    forAll(mu.boundaryField(), patchI)  // Loop over boundary patches
      {
	forAll(mu.boundaryField()[patchI], faceI)  // Loop over faces in patch
	  {
	    scalar x =  mesh.C().boundaryField()[patchI][faceI].x();
	    scalar y = mesh.C().boundaryField()[patchI][faceI].y();
	    scalar z = mesh.C().boundaryField()[patchI][faceI].z();
	    scalar a = Foam::sqrt((y * y) + (z * z)) / (rs.boundaryField()[patchI][faceI]);
	    scalar b = x / (rl.boundaryField()[patchI][faceI]);

	    
	    scalar muValue = Foam::atan2(a, b);
	    mu.boundaryFieldRef()[patchI][faceI] = muValue;

	    scalar thethaValue;

	    if (muValue <= 1e-7)
	      {
		thethaValue = 0.0;
		thetha.boundaryFieldRef()[patchI][faceI] = thethaValue;
	      }
	    else
	      {
		thethaValue = pi -  Foam::atan2(z, - y);
		thetha.boundaryFieldRef()[patchI][faceI] = thethaValue;
	      }
	  }
      }

    
    Info<< "Writing mu" << endl;
    mu.write();
    Info<< "Writing thetha" << endl;
    thetha.write();



    //initialising e_mu from the paper
    volVectorField e_mu
      (
       IOobject
       (
	"e_mu",        // Field name
	mesh.time().timeName(), // Time directory (still required but not used)
	mesh,             // Mesh reference
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
       mesh,
       vector::zero
       );
    
    //getting values of e_mu by taking the derivatice of Eq(12) w.r.t mu
    forAll(mesh.C(), cellI) {
      //NOTE I have changed the x and z values for the derivative since it
      //gave more reasonable results.

      //Below is the derivative of each component of the vector from
      //Eq (12) w.r.t mu.
      scalar e_mu_z = rs[cellI] * Foam::cos(mu[cellI]) * Foam::sin(thetha[cellI]);
      scalar e_mu_y = rs[cellI] * Foam::cos(mu[cellI]) * Foam::cos(thetha[cellI]);
      scalar e_mu_x = -rl[cellI] * Foam::sin(mu[cellI]);

      e_mu[cellI] = vector(e_mu_x, e_mu_y, e_mu_z); // Create a vector from the components
    }

    //Doing the same procedure for the boundaries
    forAll(e_mu.boundaryField(), patchI)
      {
	forAll(e_mu.boundaryField()[patchI], faceI)
	  {
	    scalar e_mu_z = rs.boundaryField()[patchI][faceI] *
	                     Foam::cos(mu.boundaryField()[patchI][faceI]) *
	                     Foam::sin(thetha.boundaryField()[patchI][faceI]);

	    scalar e_mu_y = rs.boundaryField()[patchI][faceI]
	                     *Foam::cos(mu.boundaryField()[patchI][faceI])
	                     *Foam::cos(thetha.boundaryField()[patchI][faceI]);

	    scalar e_mu_x = -rl.boundaryField()[patchI][faceI]
	                     *Foam::sin(mu.boundaryField()[patchI][faceI]);

	    // Assign the computed vector to the boundary field
	    e_mu.boundaryFieldRef()[patchI][faceI] = vector(e_mu_x, e_mu_y, e_mu_z);
	  }
      }

    e_mu /= mag(e_mu);
    Info<< "Writing e_mu" << endl;
    e_mu.write();


    //initialising e_thetha from the paper
    volVectorField e_thetha
      (
       IOobject
       (
	"e_thetha",        // Field name
	mesh.time().timeName(), // Time directory (still required but not used)
	mesh,             // Mesh reference
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
       mesh,
       vector::zero
       );
    
    //filling values of e_thetha by taking the derivative of Eq (12) w.r.t thetha
    
    forAll(mesh.C() ,cellI) {

      // The following is the derivative for the x,y,z component w.r.t thetha in Eq (12) 
      scalar e_thetha_z = rs[cellI] * Foam::sin(mu[cellI]) * Foam::cos(thetha[cellI]);
      scalar e_thetha_y = -rs[cellI] * Foam::sin(mu[cellI]) * Foam::sin(thetha[cellI]);
      scalar e_thetha_x = 0;           // 0

      e_thetha[cellI] = vector(e_thetha_x, e_thetha_y, e_thetha_z); // Create a vector from the components
    }

    //Doing the same procedure for the boundaries
    forAll(e_thetha.boundaryField(), patchI)
      {
	forAll(e_thetha.boundaryField()[patchI], faceI)
	  {
	    scalar e_thetha_z = rs.boundaryField()[patchI][faceI]
	                      * Foam::sin(mu.boundaryField()[patchI][faceI])
	                      * Foam::cos(thetha.boundaryField()[patchI][faceI]);

	    scalar e_thetha_y = -rs.boundaryField()[patchI][faceI]
	                       *Foam::sin(mu.boundaryField()[patchI][faceI])
	                       *Foam::sin(thetha.boundaryField()[patchI][faceI]);

	    scalar e_thetha_x = 0.0; // No change in the z-direction

	    // Assign the computed vector to the boundary field
	    e_thetha.boundaryFieldRef()[patchI][faceI] = vector(e_thetha_x, e_thetha_y, e_thetha_z);
	  }
      }

    e_thetha /= mag(e_thetha);
    Info<< "Writing e_thetha" << endl;
    e_thetha.write();

    
    //initialising e_t from the paper
    volVectorField e_t
      (
       IOobject
       (
	"e_t",        // Field name
	mesh.time().timeName(), // Time directory (still required but not used)
	mesh,             // Mesh reference
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
       mesh,
       vector::zero
       );

    
    //calculating the derivative of drl/dt and drs/dt
    const scalar drl_dt = r_long_epi - r_long_endo;
    const scalar drs_dt = r_short_epi - r_short_endo;

    //filling the values of e_t by taking the derivative of Eq(12) w.r.t t.
    forAll(mesh.C() ,cellI) {

      // The following is the derivative for the x,y,z component w.r.t thetha in Eq (12)
      scalar e_t_z = drs_dt * Foam::sin(mu[cellI]) * Foam::sin(thetha[cellI]);
      scalar e_t_y = drs_dt * Foam::sin(mu[cellI]) * Foam::cos(thetha[cellI]);
      scalar e_t_x = drl_dt * Foam::cos(mu[cellI]);

      e_t[cellI] = vector(e_t_x, e_t_y, e_t_z); // Create a vector from the components
    }

    //Doing the same procedure for the boundaries
    forAll(e_thetha.boundaryField(), patchI)
      {
	forAll(e_thetha.boundaryField()[patchI], faceI)
	  {
	    scalar e_t_z = drs_dt
			      * Foam::sin(mu.boundaryField()[patchI][faceI])
			      * Foam::sin(thetha.boundaryField()[patchI][faceI]);

	    scalar e_t_y = drs_dt
			       *Foam::sin(mu.boundaryField()[patchI][faceI])
			       *Foam::cos(thetha.boundaryField()[patchI][faceI]);

	    scalar e_t_x = drl_dt * Foam::cos(mu.boundaryField()[patchI][faceI]);

	    // Assign the computed vector to the boundary field
	    e_t.boundaryFieldRef()[patchI][faceI] = vector(e_t_x, e_t_y, e_t_z);
	  }
      }

    e_t /= mag(e_t);
    Info<< "Writing e_t" << endl;
    e_t.write();
    
    //Initialising the fibres
    volVectorField f1
      (
       IOobject
       (
        "f1",        // Field name
        mesh.time().timeName(), // Time directory (still required but not used)
        mesh,             // Mesh reference
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
       mesh,
       vector::zero
       );
    
    //Calculating the fibres based on Eq (15)
    forAll(mesh.C(), cellI)
    {
      f1[cellI] = 
		   ( e_mu[cellI] * Foam::sin(alphaRadians[cellI]))
		     + ( e_thetha[cellI] * Foam::cos(alphaRadians[cellI])) ;

    }

    // doing the same procedure for the boundaries
    forAll(f1.boundaryField(), patchI)
      {
	forAll(f1.boundaryField()[patchI], faceI)
	  {
	    f1.boundaryFieldRef()[patchI][faceI] = 
						    
						     (e_mu.boundaryField()[patchI][faceI]
						       *Foam::sin(alphaRadians.boundaryField()[patchI][faceI])) + 
						      (e_thetha.boundaryField()[patchI][faceI] *
						       Foam::cos(alphaRadians.boundaryField()[patchI][faceI]));
	  }
	  }

    // Optionally, write out the result to a file
    f1 /= mag(f1);
    Info<< "Writing f1" << endl;
    f1.write();

    


    volVectorField f0("f0", f1);
    f0.write();


    //Defining the sheet normal and sheet directions as per Eq (15)
    volVectorField n_("n_", e_mu ^ e_thetha );
    n_ /= mag(n_);
    Info<< "Writting n_" << endl;
    n_.write();

    volVectorField s_("s_", f0 ^ n_);
    s_ /= mag(s_);
    Info<< "Writting s_" <<endl;
    s_.write();
    
    return 0;
}


// ************************************************************************* //
