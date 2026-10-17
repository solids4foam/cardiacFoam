# leadFieldECG

The ECG of a heart+torso passive conductor from precomputed lead fields
(`buildLeadFields`): each channel is the sum, over the heart cells, of the
electrode's lead field times the heart source of the current `Vm`, plus the
constant term of the conductor's boundary data and reference.

The lead field of electrode `e` is the discrete adjoint `z_e = A^-T e_e` of
the finite-volume operator `A` of the conductor (`passiveConductor`), with
its explicit non-orthogonal, cross-diffusion and interface corrections
transposed (`leadFieldAdjoint`). The channel then equals the electrode value
of the converged forward solve `A phi = b` of the same conductor on any
mesh; `A` is not symmetric on general meshes, so the lead field of `A` itself
would not.

Electrodes lie in the torso; intracardiac electrodes (in a myocardium cell)
are not supported and stop the run with a fatal error. A bidomain myocardium
with `bathHeartPhiETrace global` is not supported either: the conductor's
heart companion is zeroGradient.
