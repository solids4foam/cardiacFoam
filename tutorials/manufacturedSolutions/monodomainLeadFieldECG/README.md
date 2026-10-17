# monodomainLeadFieldECG

Verifies `monodomainSolver` with the imported conductivity-weighted lead-field
ECG using the existing FDA manufactured voltage and an independent continuum
integral on the unit cube. The two channels have opposite polarity.

`./Allrun` generates sensitivities with `setManufacturedLeadFields` and runs
cardiacFoam. `./Allrun parallel` checks decomposed fields and collective output.
After `./Allclean`, `./Allrun fine` selects 40 cells per axis through the generated local/mesh include
for gradient/functional convergence; the default has 20 cells per axis.
Run `./regression/regressionTest.sh` with the appropriate regression build mode.

The synthetic sensitivities verify projection and source integration. They are
not clinical electrode sensitivities. `buildLeadFields` constructs those from
a separately configured passive conductor. Eikonal support is deferred.

The runner also accepts `-N`, `-dt` and `-endTime`, writing overrides to local/
without changing tracked dictionaries. Both source and observation explicitly
select the current sealed-wall gradient. Regression bounds the independent
source-field errors as well as both signed ECG functional errors.
The FDA source retains its standard diagonal propagation conductivity. For
this synthetic projection fixture, explicit physical Gi=Ge=2Gm is supplied,
so their principal-axis monodomain reduction equals the propagation tensor.
The manufactured weighted vectors remain prescribed analytical fixtures.

`-pseudoECG yes` enables a concurrent Pseudo observer while retaining ECG as
the lead observer. Its E3/E4/E5 FDA reference uses quadrature 96 versus 48 and
explicit sigmaExtracellular 0: this is an unnormalised case functional, not a
voltage-normalized lead. Regression enables both observers and checks their
independent dimensionless error cohorts and quadrature control. See
[the coexistence report](../../../docs/DUAL_ECG_VERIFICATION.md).
