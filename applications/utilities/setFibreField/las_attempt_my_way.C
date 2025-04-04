/*---------A------------------------------------------------------------------* \
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
    const scalar alphaEndo = 60.0;

    // alpha at the epicardium: SHOULD BE INPUT PARAMETER
    const scalar alphaEpi = -60.0;

    // Calculate alphaRadians

    //-A: creating a scalar field called alphaRadians. Linearly interpolating between alpha epi and endo as f varies linearly from 0 to 1 here.
    // -A: I am not sure if this is supposed to be the alpha(t) from the paper which is defined to be 90 -180 t? or if this is strictly for the rotation
    // to get the mesh.
    const volScalarField alphaRadians
    (
        "alphaRadians",
	//degToRadian(90 - 180*t)
	//(degToRad(90) - degToRad(180) * t)
        (alphaEndo*(1.0 - t) + alphaEpi*t) *(Foam::constant::mathematical::pi/180.0 )
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
     *TRYING A NEW APPROACH BY EXPLICITLY DEFINING F(U,V) BASED ON THE FUNCTION IN THE PAPER.
     *
     *-A
     */


    const volScalarField rs
      (
       "rs",
       (2.5 + 1.0*t) * 0.01 // scaling by 0.001 to get to m
       );
    Info<< "Writing rs" << endl;
    rs.write();

    const volScalarField rl
      (
       "rl",
       (9.0 + 0.7*t) * 0.01 // scaling by 0.001 to get in m
       );
    Info<< "Writing rl" << endl;
    rl.write();


    //finding the minimum z-value to change the computation close to the apex
    scalar minZ = 100;  // Initialize with a large value

    forAll(mesh.C(), cellI)
      {
	scalar z = mesh.C()[cellI].z();
	minZ = min(minZ,z);

      };

    /*
    scalar x = mesh.C().x();
    scalar y = mesh.C().y();
    scalar z = mesh.C().z();

    scalar a = Foam::sqrt((y * y ) + (z * z ) ) / rs;
    scalar b = x / rl;
    scalar u =*/
    //const vectorField& C = mesh.C();

    //scalarField x("x",C.X());
    //scalarField y("y", C.Y());
    //scalarField z("z", C.Z());



    // Access the mesh (directly use 'mesh' if you are inside a solver)
    // Get the mesh points (coordinates)
    const pointField& points = mesh.points();

    // Initialize fields to hold x, y, and z values
    volScalarField x
      (
       IOobject
       (
        "x",  // Field name
        runTime.timeName(),  // Time step
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
       mesh,
       dimensionedScalar("x", dimLength, 0.0) // Initialize with 0
       );

    volScalarField y
      (
       IOobject
       (
        "y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
       mesh,
       dimensionedScalar("y", dimLength, 0.0) // Initialize with 0
       );

    volScalarField z
      (
       IOobject
       (
        "z",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
       mesh,
       dimensionedScalar("z", dimLength, 0.0) // Initialize with 0
       );

    // Loop through all the points and assign x, y, z values
    for (label i = 0; i < points.size(); ++i)
      {
	x[i] = points[i].x();  // X component
	y[i] = points[i].y();  // Y component
	z[i] = points[i].z();  // Z component
      }

    // Now you can use x, y, z for further calculations or processing

    
    
    volScalarField a("a", Foam::sqrt((y * y) + (z * z)) / rs);
    volScalarField b("b", x / rl);
    volScalarField u("u", Foam::atan2(a, b));
    volScalarField v("v", Foam::constant::mathematical::pi - Foam::atan2(z, -y));
    v = pos(u - scalar(1e-7)) * v + neg(u - scalar(1e-7)) * scalar(0);

    Info<< "Writing u" << endl;
    u.write();
    Info<< "Writing v" << endl;
    v.write();
    /*volVectorField dx_du(
			 "dx_du",

			 vectorField
			 (

			  -rl * Foam::sin(u),

			  rs * Foam::cos(u) *Foam::cos(v),

			  rs * Foam::cos(u) * Foam::sin(v)
			  )
			 );
    dx_du /= mag(dx_du);
    */

    const volScalarField dxu
      (
       "dxu",
       -(rl * Foam::sin(u)) // scaling by 0.001 to get in m                                                                                                 
       );
    
    Info<< "Writing dxu" << endl;
    dxu.write();

    const volScalarField dyu
      (
       "dyu",
       (rs * Foam::cos(u) * Foam::cos(v)) // scaling by 0.001 to get in m
       );

    Info<< "Writing dyu" << endl;
    dyu.write();

    const volScalarField dzu
      (
       "dzu",
       (rs * Foam::cos(u) * Foam::sin(v))
       );

    Info<< "Writting dzu" << endl;
    dzu.write();
    

    volVectorField dx_du("dx_du", et);

    forAll(et, cellI) {
      // Calculate the components of d𝐱/dv for each i-th sample
      scalar dx_du_x = dxu[cellI]; // 𝑟𝑠 * sin(u) * cos(v)
      scalar dx_du_y = dyu[cellI];
      scalar dx_du_z = dzu[cellI];

      // Assign the value to the dx_dv field (as a vector)
      dx_du[cellI] = vector(dx_du_x, dx_du_y, dx_du_z); // Create a vector from the components
    }
    dx_du /= mag(dx_du);

    Info<< "Writting dx_du" << endl;
    dx_du.write();

    const volScalarField dxv
      (
       "dxv",
       -(rl * Foam::sin(u)) * 0.0         

       );

    Info<< "Writing dxv" << endl;
    dxv.write();

    const volScalarField dyv
      (
       "dyv",
       (-rs * Foam::sin(u) * Foam::sin(v)) // scaling by 0.001 to get in m
       );

    Info<< "Writing dyv" << endl;
    dyv.write();

    const volScalarField dzv
      (
       "dzv",
       (rs * Foam::sin(u) * Foam::cos(v))
       );

    Info<< "Writting dzv" << endl;
    dzv.write();

    volVectorField dx_dv("dx_dv", et);

    forAll(et, cellI) {
      // Calculate the components of d𝐱/dv for each i-th sample
      scalar dx_dv_x = dxv[cellI]; // 𝑟𝑠 * sin(u) * cos(v)
      scalar dx_dv_y = dyv[cellI];
      scalar dx_dv_z = dzv[cellI];

      // Assign the value to the dx_dv field (as a vector)
      dx_dv[cellI] = vector(dx_dv_x, dx_dv_y, dx_dv_z); // Create a vector from the components
    }
    dx_dv /= mag(dx_dv);

    Info<< "Writting dx_dv" << endl;
    dx_dv.write();
    
    /*volVectorField dx_du
      (
       IOobject
       (
        "dx_du",                      // Field name
        runTime.timeName(),            // Time step
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
	),
       mesh,
       vectorField
       (
        dxu,  // X component
        dyu,  // Y component
        dzu   // Z component
	)
       );*/

    // Normalize the vector field
    volVectorField f1(
		      "f1",
		      Foam::sin(alphaRadians) * dx_du + Foam::cos(alphaRadians) * dx_dv
		      );
    f1 /= mag(f1);
    Info<< "Writing f1" << endl;
    f1.write();

    volVectorField ns("ns", dx_du ^ dx_dv );
    ns /= mag(ns);
    Info<< "Writing ns" << endl;
    ns.write();

    volVectorField s_d("s_d", f1 ^ ns);
    s_d /= mag(s_d);

    Info<< "Writing s_d" << endl;
    s_d.write();
    
    /*
    scalar threshold = minZ;
    volVectorField uv("uv", el);
    forAll(et, cellI)
      {
	scalar x = mesh.C()[cellI].x();
	scalar y = mesh.C()[cellI].y();
	scalar z = mesh.C()[cellI].z();

	//initialising u and v.
	scalar uValue;
	scalar vValue;
	if(z<=threshold + 0.000001 )
	  {
	    scalar a = (Foam::sqrt((y * y) + (z * z))) /(rs[cellI]);
	    scalar b = x/rl[cellI];
	    uValue = Foam::atan(a/b);
	    uv[cellI].x() = uValue;
	    uv[cellI].z()=0;
	    //changing v-values based on u values threshold.
	    if(uValue <= 1e-7)
	      {
		vValue = 0;
		uv[cellI].y() = vValue;
	      }
	    else
	      {
		vValue = Foam::constant::mathematical::pi - Foam::atan(-z/y);
		uv[cellI].y() = vValue;
	      }


	  }
	else
	  {
	    uValue = Foam::acos(x/ rl[cellI]);
	    uv[cellI].x() = uValue;

	    vValue = Foam::acos(y / (rs[cellI] * Foam::sin(uValue)));
	    uv[cellI].y() = vValue;
	    uv[cellI].z() =0;
	  }

      }
    
    forAll(et.boundaryField(), patchI) // same idea as above but applied to the boundary.                                                                  
    {
      const word patchName = mesh.boundaryMesh()[patchI].name();
      forAll(et.boundaryField()[patchI], cellI)
	{
	  scalar x = mesh.Cf()[cellI].x();
	  scalar y = mesh.Cf()[cellI].y();
	  scalar z = mesh.Cf()[cellI].z();
	  scalar umin= -Foam::constant::mathematical::pi, umax= Foam::constant::mathematical::pi;
	  scalar vmin = -Foam::constant::mathematical::pi, vmax = Foam::constant::mathematical::pi;

	  if(patchName == "inside")
	    {
	      umin = -Foam::constant::mathematical::pi;
	      umax = -Foam::acos(5.0/17.0);
	      vmin = -Foam::constant::mathematical::pi;
	      vmax = Foam::constant::mathematical::pi;
	      
	    }
	  
	  else if(patchName == "outside")
	    {
	      umin = -Foam::constant::mathematical::pi;
              umax = -Foam::acos(5.0/20.0);
              vmin = -Foam::constant::mathematical::pi;
              vmax = Foam::constant::mathematical::pi;
	    }

	  else
	    {
	      continue;
	    }

	  
	  //initialising u and v.
	  scalar uValue;
	  scalar vValue;
	  if(z<=threshold + 0.00000001 )
	    {
	      scalar a = (Foam::sqrt((y * y) + (z * z))) /(rs[cellI]);
	      scalar b = x/rl[cellI];
	      uValue = Foam::clamp(Foam::atan(a/b), umin, umax);
	      uv.boundaryFieldRef()[patchI][cellI].x() = uValue;
	      uv.boundaryFieldRef()[patchI][cellI].z()=0;

	      //changing v-values based on u values threshold.
	      if(uValue <= 1e-7)
		{
		  vValue = 0;
		  uv.boundaryFieldRef()[patchI][cellI].y() = vValue;
		}
	      else
		{
		  vValue = Foam::clamp(Foam::constant::mathematical::pi - Foam::atan(-z/y), vmin, vmax);
		  uv.boundaryFieldRef()[patchI][cellI].y() = vValue;
		}


	    }
	  else
	    {
	      uValue = Foam::clamp(Foam::acos(x/ rl[cellI]), umin, umax);
	      uv.boundaryFieldRef()[patchI][cellI].x() = uValue;

	      vValue = Foam::clamp(Foam::acos(y / (rs[cellI] * Foam::sin(uValue))), vmin, vmax);
	      uv.boundaryFieldRef()[patchI][cellI].y() = vValue;
	      uv.boundaryFieldRef()[patchI][cellI].z() =0;
	    }
	}
    }
    
    Info << "Writing uv" << endl;
    uv.write(); */

    /*
    volVectorField dx_du
      (
       IOobject
       (
	"dx_du",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       mesh,
       dimensionedVector("zero", dimLength, vector(0, 0, 0))
       );
    */

    /*
    volVectorField dx_du("dx_du", et);

    forAll(et, cellI) {
      // Calculate the components of d𝐱/dv for each i-th sample
      scalar dx_du_z = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::cos(uv[cellI].y()); // 𝑟𝑠 * sin(u) * cos(v)
      scalar dx_du_y = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::sin(uv[cellI].y()); // 𝑟𝑠 * sin(u) * sin(v)
      scalar dx_du_x = rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::cos(uv[cellI].y());  // 0

      // Assign the value to the dx_dv field (as a vector)
      dx_du[cellI] = vector(dx_du_x, dx_du_y, dx_du_z); // Create a vector from the components
    }
    forAll(el.boundaryField(), patchI)
      {
	forAll(el.boundaryField()[patchI], cellI)
	  {
	    scalar dx_du_z = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::cos(uv[cellI].y()); // 𝑟𝑠 * sin(u) * cos(v)
	    scalar dx_du_y = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::sin(uv[cellI].y()); // 𝑟𝑠 * sin(u) * sin(v)
	    scalar dx_du_x = rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::cos(uv[cellI].y());  // 0
	    dx_du.boundaryFieldRef()[patchI][cellI] = vector(dx_du_x, dx_du_y, dx_du_z);
	  }
      }
    dx_du = dx_du/mag(dx_du);
     
    Info << "Writing dx_du field with updated boundary conditions" << endl;
    dx_du.write();

    */
    /*
    volVectorField dx_dv
      (
       IOobject
       (
        "dx_dv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedVector("zero", dimLength, vector(0, 0, 0))
       );

    */
    /*
    volVectorField dx_dv("dx_dv", et);
    forAll(et ,cellI) {
      // Calculate the components of d𝐱/du for each i-th sample                                                                                 
      scalar dx_dv_z = -rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::sin(uv[cellI].y()); // 𝑟𝑠 * cos(u) * cos(v)  
      scalar dx_dv_y = rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::cos(uv[cellI].y()); // 𝑟𝑠 * cos(u) * sin(v)
      scalar dx_dv_x = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::sin(uv[cellI].y());
      dx_dv[cellI] = vector(dx_dv_x, dx_dv_y, dx_dv_z); // Create a vector from the components
    }

    
    forAll(el.boundaryField(), patchI)
      {
	forAll(el.boundaryField()[patchI], cellI)
	  {
	    scalar dx_dv_z = -rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::sin(uv[cellI].y()); // 𝑟𝑠 * cos(u) * cos(v)
	    scalar dx_dv_y = rs[cellI] * Foam::sin(uv[cellI].x()) * Foam::cos(uv[cellI].y()); // 𝑟𝑠 * cos(u) * sin(v)
	    scalar dx_dv_x = rs[cellI] * Foam::cos(uv[cellI].x()) * Foam::sin(uv[cellI].y());
	    dx_dv.boundaryFieldRef()[patchI][cellI] = vector(dx_dv_x, dx_dv_y, dx_dv_z);
	    
	  }
      }
    dx_dv = dx_dv/mag(dx_dv);
    
    Info<< "Writing dx_dv" << endl;
    dx_dv.write();

    */

    /*
    volVectorField f1
      (
       IOobject
       (
	"f1",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       mesh,
       dimensionedVector("zero", dimLength, vector(0, 0, 0))
       );
    */

    //volVectorField el("el", en ^ et);
    //el /= mag(el);
    
    //volVectorField f1("f1", et);
    /*
    forAll(et, cellI)
      {
	scalar norm_dx_du = mag(dx_du[cellI]);
	scalar norm_dx_dv = mag(dx_dv[cellI]);
	f1[cellI] = ( Foam::sin(alphaRadians[cellI]) * (dx_du[cellI]) + ( Foam::cos(alphaRadians[cellI]) * (dx_dv[cellI])));

      }
    forAll(et.boundaryField(), patchI)
      {
	forAll(et.boundaryField()[patchI], cellI){
	  scalar norm_dx_du = mag(dx_du[cellI]);
	  scalar norm_dx_dv = mag(dx_dv[cellI]);
	  f1.boundaryFieldRef()[patchI][cellI] = ( Foam::sin(alphaRadians[cellI]) * (dx_du[cellI])) + ( Foam::cos(alphaRadians[cellI]) * (dx_dv[cellI]));
	  //scalar norm_f1 = mag(f1[cellI]);
	  //f1.boundaryFieldRef()[patchI][cellI] = f1[cellI] / norm_f1;
	    
	}
	}*/
    //volVectorField f1 = Foam::sin(alphaRadians) * dx_du + Foam::cos(alphaRadians) * dx_dv;

    /*
    volVectorField f1( Foam::sin(alphaRadians) * dx_du + Foam::cos(alphaRadians) * dx_dv );
    f1 /= mag(f1);

    // Optionally, write out the result to a file
    Info<< "Writing f1" << endl;
    f1.write();

    volVectorField ns("ns", dx_du ^ dx_dv);                                                                                                                     
    ns /= mag(ns);
    Info<< "Writing ns" << endl;
    ns.write();

    volVectorField s_d("s_d", f1 ^ ns);
    s_d /= mag(s_d);
    Info<< "Writing s_d" << endl;
    s_d.write();

    */


    
    /*End of trying new approach */
    
    //initialising the f0 to el field
    volVectorField f0("f0", el);
    //looping though the values of et for each cell.
    forAll(et, cellI)
    {
        const tensor rotT = Ra(et[cellI], alphaRadians[cellI]); //computing the rotation matrix for.. rotation around the longitudinal by alpharadians. 
        f0[cellI] = transform(-rotT, f0[cellI]); // Applying the rotation.
    }
    /*forAll(et.boundaryField(), patchI) // same idea as above but applied to the boundary.
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
	}*/

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
