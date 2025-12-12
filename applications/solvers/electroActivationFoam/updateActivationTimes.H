// Update activationTimes field (used in Niederer benchmark)
const scalarField& VmI = Vm.primitiveFieldRef();
const scalarField& VmOldI = Vm.oldTime().primitiveFieldRef();
scalarField& activationTimeI = activationTime.primitiveFieldRef();
const scalar oldTime = runTime.value() - runTime.deltaTValue();
forAll(activationTimeI, cellI)
{
    if (calculateActivationTime[cellI])
    {
        if (VmI[cellI] > SMALL)
        {
            calculateActivationTime[cellI] = false;

            // Use current time
            // activationTimeI[cellI] = runTime.value();

            // Linearly interpolate for more accuracy
            const scalar w = (0.0 - VmOldI[cellI])/(VmI[cellI] - VmOldI[cellI]);

            activationTimeI[cellI] = oldTime + w*runTime.deltaTValue();
        }
    }
}
activationTime.correctBoundaryConditions();
