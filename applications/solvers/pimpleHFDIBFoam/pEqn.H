volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Immersed flux
// Philip: should we instead include this in phiHbyA?
const volVectorField FbyA(rAU*f);
const surfaceScalarField phiFbyA("phiFbyA", fvc::flux(FbyA));

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) + fvc::div(phiFbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(p.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p) + rAU*f;
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

// Update immersed fields
HFDIBDEM.updateFSCoupling(lambda,f);

Ui.correctBoundaryConditions();
lambda.correctBoundaryConditions();
f.correctBoundaryConditions();
refineF.correctBoundaryConditions();
surface.correctBoundaryConditions();

HFDIBDEM.interpolateIB(U,Ui,lambda);

Ui.correctBoundaryConditions();
lambda.correctBoundaryConditions();
f.correctBoundaryConditions();
refineF.correctBoundaryConditions();
surface.correctBoundaryConditions();

forAll(surface, sI)
{
    if (lambda[sI] > SMALL) //0)
    {
        surface[sI] = 1;
    }
    else
    {
        surface[sI] = 0;
    }
}
Ui.correctBoundaryConditions();
lambda.correctBoundaryConditions();
f.correctBoundaryConditions();
refineF.correctBoundaryConditions();
surface.correctBoundaryConditions();

// Update immersed force to penalise the difference between U and Ui
f += ibCouplingCoeff*surface*(Ui - U)/runTime.deltaT();
