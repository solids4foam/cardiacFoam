# leadFieldECG

The lead-field evaluation of `torsoECG` with `torsoCoupling leadField`: the
ECG of a heart+torso passive conductor from precomputed lead fields
(`buildLeadFields`). Each channel is the sum, over the heart cells, of the
electrode's lead field times the heart source of the current `Vm`, plus the
constant term of the conductor's boundary data and reference.

```
torso
{
    cellZones         ( bath );
    conductivityField bodyAndOrgansConductivity;
    interfaceConductivityInterpolation distanceWeightedHarmonic;
    groundPatches     { ... }     // or phiERefPoint and phiEReferenceValue
}

torsoCoupling leadField;

ecgDomains
{
    ECG
    {
        ecgSolver          torsoECG;
        electrodePositions { ... }
        // monodomainSolver only: Gi and Ge with Gi (Gi+Ge)^-1 Ge equal to
        // its conductivity
        conductivitySource uniform;
        conductivityIntracellular [ -1 -3 3 0 0 2 0 ] ( ... );
        conductivityExtracellular [ -1 -3 3 0 0 2 0 ] ( ... );
    }
}
```

`leadFieldConductor` builds the conductor from these blocks, `leadFieldAdjoint`
transposes its explicit part for `buildLeadFields`, and `leadFieldECG` loads
and evaluates the lead fields.

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
