# checkLeadFieldRecovery

Controlled verification of reciprocal vector projection against production
passive extracellular recovery and `torsoECG` cell observations. Requires the
orthogonal embedded-heart column used in configurationChecks, physical Gi/Ge,
and precomputed sensitivities. This utility is not a generic anatomical
recovery tool: its absolute 1e-8 V criterion is specific to this control.

It constructs a verification-owned bidomain domain, prescribes a nonzero
linear Vm snapshot, and iterates the production `solvePhiE()` with Vm held
fixed. It combines sampled torso potentials with the exact balanced weights
from `leadFieldProperties` and compares against projection of the same Vm.
Both Vm internal and boundary arrays must remain unchanged. Manufactured
source forcing is forbidden. No ionic step or physical time advance occurs.

Serial and MPI configuration hooks exercise production ownership and solve
paths. Agreement on this orthogonal scalar problem does not establish an
exact adjoint of the corrected, anisotropic or nonorthogonal operator.
