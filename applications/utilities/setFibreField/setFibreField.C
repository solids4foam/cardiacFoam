/*--------A---A------A---------------------------------------------------------*\
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
    #include "vector.H"
    #include "tensor.H"
    #include "cmath"

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
    const scalar fibre_Angle_Endo = 90.0;

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar fibre_Angle_Epi = -90.0;

    const scalar pi = Foam::constant::mathematical::pi;

    const volScalarField alphaRadians
    (
        "alphaRadians",
        (fibre_Angle_Endo + (fibre_Angle_Epi - fibre_Angle_Endo)*t) * pi/180
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
    //from Eq (16) of the paper we get rs and rl
    const volScalarField rs
    (
        "rs",
	7 + 3*t
    );
    Info<< "Writing rs" << endl;
    rs.write();

    
    const volScalarField rl
    (
        "rl",
        17 + 3*t
    );
    Info<< "Writing rl" << endl;
    rl.write();

    
    //initialising uu from the paper
    volScalarField uu
      (
       IOobject
       (
        "uu",                    // Name of the field
        runTime.timeName(),     // Time directory (e.g., "0", "1", etc.)
        mesh,                   // Mesh to associate the field with
        IOobject::NO_READ,      // Do not read from disk (no existing file)
        IOobject::AUTO_WRITE    // Write the field automatically to disk
	),
       mesh,
       dimensionedScalar("uu", dimless, 0.0) // Mesh to define the field on
       );

    //Initialising vv from the paper
    volScalarField vv
      (
       IOobject
       (
        "vv",                    // Name of the field                                                                                                        
        runTime.timeName(),     // Time directory (e.g., "0", "1", etc.)                                                                                    
        mesh,                   // Mesh to associate the field with                                                                                         
        IOobject::NO_READ,      // Do not read from disk (no existing file)                                                                                 
        IOobject::AUTO_WRITE    // Write the field automatically to disk                                                                                    
        ),
       mesh,
       dimensionedScalar("vv", dimless, 0.0)                                                                                             
       );

    volScalarField q
      (
       IOobject
       (
	"q",

	runTime.timeName(),

	mesh,

	IOobject::NO_READ,

	IOobject::AUTO_WRITE

	),
       mesh,
       dimensionedScalar("q", dimless, 0.0)
       );

    //looping through an updating the mu and thetha values based procedure described in paper
    forAll(mesh.C(), cellI)
      {
        
	scalar x =  mesh.C()[cellI].x();
	scalar y = mesh.C()[cellI].y();
	scalar z = mesh.C()[cellI].z();

	vv[cellI] = Foam::atan2(-y, -x);
	
	if(std::abs(Foam::cos(vv[cellI])) > 1e-6)
	  {
	    q[cellI] = x / Foam::cos(vv[cellI]);
	  }
	
	else
	  {
	    q[cellI] = y / Foam::sin(vv[cellI]);
	  }

	uu[cellI] = Foam::acos(z / rl[cellI]);

	if(uu[cellI] > 0)
	  {
	    uu[cellI] = -uu[cellI];
	  }
      }

    forAll(mesh.C().boundaryField(), patchI)  // Loop over boundary patches
      {
	forAll(mesh.C().boundaryField()[patchI], faceI)  // Loop over faces in patch
	  {
	    scalar x =  mesh.C().boundaryField()[patchI][faceI].x();
	    scalar y = mesh.C().boundaryField()[patchI][faceI].y();
	    scalar z = mesh.C().boundaryField()[patchI][faceI].z();

	    vv.boundaryFieldRef()[patchI][faceI] = Foam::atan2(-y, -x);
	    
	    if (std::abs(Foam::cos(vv.boundaryFieldRef()[patchI][faceI])) > 1e-6)
	      {
		q.boundaryFieldRef()[patchI][faceI] = x / Foam::cos(vv.boundaryField()[patchI][faceI]);
	      }
	    else
	      {
		q.boundaryFieldRef()[patchI][faceI] = y / Foam::sin(vv.boundaryField()[patchI][faceI]);
	      }


	    uu.boundaryFieldRef()[patchI][faceI] = Foam::acos(z / rl.boundaryField()[patchI][faceI]);

	    if(uu.boundaryFieldRef()[patchI][faceI] > 0)
	      {
		uu.boundaryFieldRef()[patchI][faceI] = -uu.boundaryField()[patchI][faceI];
	      }
	  }
      }

    
    Info<< "Writing uu" << endl;
    uu.write();
    Info<< "Writing vv" << endl;
    vv.write();
    Info<< "Writing q" << endl;
    q.write();
    
    
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
      scalar u = uu[cellI];
      scalar v = vv[cellI];
      scalar fibre_angle = alphaRadians[cellI];
      scalar short_r = rs[cellI];
      scalar long_r = rl[cellI];


      vector deriv_dir(Foam::sin(fibre_angle), Foam::cos(fibre_angle), 0.0);


      vector Mcol0(
		   short_r * std::cos(u) * std::cos(v),
		   short_r * std::cos(u) * std::sin(v),
		   -long_r * std::sin(u)
		   );
      vector Mcol1(
		   -short_r * std::sin(u) * std::sin(v),
		   short_r * std::sin(u) * std::cos(v),
		   0
		   );

      Mcol0 /= mag(Mcol0);
      Mcol1 /= mag(Mcol1);

      tensor M(
	       Mcol0.x(), Mcol1.x(), 0,
	       Mcol0.y(), Mcol1.y(), 0,
	       Mcol0.z(), Mcol1.z(), 0
	       );

      vector fVal = M & deriv_dir;

      f1[cellI] = fVal;
    }

    // doing the same procedure for the boundaries
    forAll(f1.boundaryField(), patchI)
      {
	forAll(f1.boundaryField()[patchI], faceI)
	  {
	    scalar u = uu.boundaryField()[patchI][faceI];
	    scalar v = vv.boundaryField()[patchI][faceI];
	    scalar fibre_angle = alphaRadians.boundaryField()[patchI][faceI];
	    scalar short_r = rs.boundaryField()[patchI][faceI];
	    scalar long_r = rl.boundaryField()[patchI][faceI];
	    


	    vector deriv_dir(Foam::sin(fibre_angle), Foam::cos(fibre_angle), 0.0);


	    vector Mcol0(
			 short_r * std::cos(u) * std::cos(v),
			 short_r * std::cos(u) * std::sin(v),
			 -long_r * std::sin(u)
			 );
	    vector Mcol1(
			 -short_r * std::sin(u) * std::sin(v),
			 short_r * std::sin(u) * std::cos(v),
			 0
			 );

	    Mcol0 /= mag(Mcol0);
	    Mcol1 /= mag(Mcol1);

	    tensor M(
		     Mcol0.x(), Mcol1.x(), 0,
		     Mcol0.y(), Mcol1.y(), 0,
		     Mcol0.z(), Mcol1.z(), 0
		     );

	    vector fVal = M & deriv_dir;

	    f1.boundaryFieldRef()[patchI][faceI] = fVal;
	    
	    
	  }
	  }

    // Optionally, write out the result to a file
    f1 /= mag(f1);
    Info<< "Writing f1" << endl;
    f1.write();

    


    volVectorField f0("f0", f1);
    f0.write();

    
    
    return 0;
}


// ************************************************************************* //
