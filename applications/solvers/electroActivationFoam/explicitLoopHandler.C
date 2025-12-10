#include "explicitLoopHandler.H"
#include "tmanufacturedFDA.H"

explicitLoopHandler::explicitLoopHandler
(
    const fvMesh& mesh,
    ionicModel& model
)
:
    // mesh_(mesh),
    ionicModel_(model),
    dx_(0.0),
    dim_(0),
    Def_max_(0.0)
{}



//----------------------------------------------------------------------
// Compute stability limit dt <= CFL * dx^2 / (dim * Def_max)
//----------------------------------------------------------------------
scalar explicitLoopHandler::computeStableDt(const scalar CFL) const
{
    if (dim_ == 0 || Def_max_ <= SMALL)
        return GREAT;

    return CFL * dx_ * dx_ / (dim_ * Def_max_);
}


//----------------------------------------------------------------------
// Initialization: compute Def_max, stable dt
//----------------------------------------------------------------------
void explicitLoopHandler::initializeExplicit
(
    scalar& dt,
    int& nsteps,
    const scalar chiVal,
    const scalar CmVal,
    const volTensorField& conductivity,
    const scalar CFL,
    const scalar dx,
    const int dim
)
{
    // Store caller-provided dx and dim (includes MS corrections)
    dx_  = dx;
    dim_ = dim;

    // Extract max diffusion coefficient
    scalar Dxx = max(conductivity.component(tensor::XX)).value();
    scalar Dyy = max(conductivity.component(tensor::YY)).value();
    scalar Dzz = max(conductivity.component(tensor::ZZ)).value();

    scalar Dmax_raw = max(Dxx, max(Dyy, Dzz));

    // Effective diffusion coefficient: D / (chi * Cm)
    Def_max_ = Dmax_raw / (chiVal * CmVal);

    // Compute stable dt
    scalar dt_stable = computeStableDt(CFL);

    dt = min(dt, dt_stable);
    reduce(dt, maxOp<scalar>());

    Info << "ExplicitLoop: dx = " << dx_
         << ", dim = " << dim_
         << ", Def_max = " << Def_max_
         << ", stable dt = " << dt_stable
         << ", using dt = " << dt << nl;

    // nsteps recomputed in solver; we leave it untouched here
    (void)nsteps;
}

//----------------------------------------------------------------------
// Perform a single explicit step
//----------------------------------------------------------------------
void explicitLoopHandler::explicitLoop
(
    const scalar t0,
    const scalar dt,
    volScalarField& Vm,
    volScalarField& Iion,
    Field<Field<scalar>>& states,
    volScalarField& externalStimulusCurrent,
    const labelList& stimulusCellIDs,
    const scalar stimulusIntensity,
    const scalar stimulusDuration,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const volTensorField& conductivity
)
{
    // // 0) Manufactured-specific: store "old" internal states
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).updateStatesOld();
    }

    // 1) External stimulus (internal field only)
    scalarField& externalStimulusCurrentI = externalStimulusCurrent;
    externalStimulusCurrentI = 0.0;

    if (t0 <= stimulusDuration)
    {
        forAll(stimulusCellIDs, cI)
        {
            const label id = stimulusCellIDs[cI];
            externalStimulusCurrentI[id] = stimulusIntensity;
        }
    }
    externalStimulusCurrent.correctBoundaryConditions();

    // 2) Ionic current using OLD Vm
    ionicModel_.calculateCurrent
    (
        t0,
        dt,
        Vm.internalField(),   // Vm_old
        Iion,
        states
    );
    Iion.correctBoundaryConditions();

    if (ionicModel_.hasManufacturedSolution())
    {
        Iion /= Cm.value();
    }

    solve
    (
        chi*Cm * fvm::ddt(Vm)
        == fvm::laplacian(conductivity, Vm)
        - chi*Cm*Iion
        + externalStimulusCurrent
    );

    Vm.correctBoundaryConditions();

//4) Manufactured-specific: reset internal states back to OLD
    if (ionicModel_.hasManufacturedSolution())
    {
        refCast<tmanufacturedFDA>(ionicModel_).resetStatesToStatesOld();
    }

// 5) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
    ionicModel_.solveODE
    (
        t0,
        dt,
        Vm.internalField(),    // Vm_new
        Iion,
        states
    );

}


