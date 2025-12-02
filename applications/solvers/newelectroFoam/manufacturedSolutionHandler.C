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
    ionicModel_(model),   // <--- reference bound here
    dx_(0.0),
    dim_(0),
    N_(1)
{}


// 1. Manufactured-solution INITIALIZATION
// ===============================================================
void manufacturedSolutionHandler::initializeManufactured
(
    volScalarField& Vm,
    volScalarField& u1,
    volScalarField& u2,
    volScalarField& u3
)
{
    Info << "Initializing manufactured-solution analytic fields..." << nl;

    // Analytic initialization from tmanufacturedFields
    ionicModel_.initializeFields(Vm, u1, u2, u3, mesh_.C());

    // Manufactured-solution domain dimension
    dim_ = ionicModel_.tissue();

    const double totalCells =
        returnReduce(mesh_.nCells(), sumOp<int>());

    N_ = std::round(Foam::pow(totalCells, 1.0/dim_));

    dx_ = 1.0 / scalar(N_);

    Info << "MS initialized: totalCells=" << totalCells
         << " N=" << N_
         << " dx=" << dx_
         << " dim=" << dim_
         << nl;
}


// 2. NON-manufactured INITIALIZATION
// ===============================================================
void manufacturedSolutionHandler::initializeNonManufactured()
{
    dx_ = Foam::cbrt(mesh_.V().average().value());
    dim_ = mesh_.nGeometricD();

    Info << "Non-MS run: structured avg dx = " << dx_
         << ", dim = " << dim_ << nl;
}


// 3. MS POST-PROCESSING
// ===============================================================
void manufacturedSolutionHandler::postProcess
(
    const volScalarField& Vm,
    const volScalarField& u1,
    const volScalarField& u2,
    const volScalarField& u3,
    const scalar dt,
    const int nsteps,
    const bool solveExplicit
) const
{
    Info << "\nCalculating manufactured-solution errors..." << nl;

    scalarField x = mesh_.C().component(vector::X);
    scalarField y = mesh_.C().component(vector::Y);
    scalarField z = mesh_.C().component(vector::Z);

    volScalarField VmMS(Vm); VmMS.rename("VmMS");
    volScalarField u1MS(u1); u1MS.rename("u1MS");
    volScalarField u2MS(u2); u2MS.rename("u2MS");
    volScalarField u3MS(u3); u3MS.rename("u3MS");

    // Manufactured-solution export (analytic fields)
    refCast<tmanufacturedFDA>(ionicModel_).exportManufacturedStates
    (
        VmMS, u1MS, u2MS, u3MS
    );

    const scalar Tfinal = Vm.time().value();

    computeAndPrintErrors
    (
        VmMS.internalField(),
        u1MS.internalField(),
        u2MS.internalField(),
        x, y, z,
        Tfinal,
        ionicModel_.tissue(),
        N_, dx_, dt, nsteps,
        solveExplicit
    );
}
