#include "manufacturedSolutionHandler.H"
#include "tmanufacturedFields.H"


// Constructor
manufacturedSolutionHandler::manufacturedSolutionHandler
(
    const fvMesh& mesh,
    ionicModel& model     
)
:
    mesh_(mesh),
    ionicModel_(model),  
    dx_(0.0),
    dim_(0),
    N_(1)
{}


// 1. Manufactured-solution INITIALIZATION
// ===============================================================
void manufacturedSolutionHandler::initializeManufactured
(
    volScalarField& Vm,
    List<volScalarField*>& outFields,
    scalar dxstructured,
    int dim
)
{
    Info << "Manufactured solution needs to export u1,u2,u3 always for error computation" << nl;

    // Get exported variable names from the ionic model
    Foam::wordList names = ionicModel_.exportedFieldNames();
    
    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");
    const label iu3 = names.find("u3");

    volScalarField& u1 = *outFields[iu1];
    volScalarField& u2 = *outFields[iu2];
    volScalarField& u3 = *outFields[iu3];

    // Analytic initialization from tmanufacturedFields via the ionic model
    ionicModel_.initializeFields(Vm, u1, u2, u3, mesh_.C());
    dim_ = dim;     // from the mesh-based general dimension
    dx_  = dxstructured;      // the generalized dx you computed in main()      
    const double totalCells =
         returnReduce(mesh_.nCells(), sumOp<int>());

    N_ = std::round(Foam::pow(totalCells, 1.0/dim_));
}


// 3. MS POST-PROCESSING
// ===============================================================
void manufacturedSolutionHandler::postProcess
(
    const volScalarField& Vm,
    const List<volScalarField*>& outFields,
    const scalar dt,
    const int nsteps,
    const bool solveExplicit
) const
{
    Info << "\nCalculating manufactured-solution errors..." << nl;

    // Lookup u1 and u2 from exported fields (numerical solution)
    Foam::wordList names = ionicModel_.exportedFieldNames();

    const label iu1 = names.find("u1");
    const label iu2 = names.find("u2");
    
    const volScalarField& u1 = *outFields[iu1];
    const volScalarField& u2 = *outFields[iu2];

    scalarField x = mesh_.C().component(vector::X);
    scalarField y = mesh_.C().component(vector::Y);
    scalarField z = mesh_.C().component(vector::Z);

    const scalar Tfinal = Vm.time().value();

    // Compare numerical Vm,u1,u2 with analytic manufactured solution
    computeAndPrintErrors
    (
        Vm.internalField(),
        u1.internalField(),
        u2.internalField(),
        x, y, z,
        Tfinal,
        ionicModel_.tissue(),
        N_, dx_, dt, nsteps,
        solveExplicit
    );
}
