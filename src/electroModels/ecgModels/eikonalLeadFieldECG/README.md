# eikonalLeadFieldECG

Passive reciprocal ECG recording for `eikonalSolver`, using a spatially
uniform action-potential template and precomputed lead sensitivities.

```foam
ecgDomains
{
    Reciprocal
    {
        ecgSolver eikonalLeadFieldECG;
        templateType endo;
        sampling { start 0; end 1; deltaT 0.001; }
    }
}
```

The source is `grad(Vm) = -U'(t-tau) grad(tau)`. The observer uses a fresh
private activation-time gradient, the source owner's wall policy, and the
shared template helper. The `eikonalLeadFieldGradTau` entry selects its
gradient scheme, falling back to the case default.

The existing `leadFieldData` importer loads the sensitivities and validates
mesh/provenance and physical Gi/Ge. Weighted sensitivities already contain
intracellular conductivity and normalized reciprocal current. Output is
`postProcessing/eikonalLeadFieldECG.dat`, in seconds and volts, sampled once
per observer instance independently of propagation steps.

Supported templates are `endo`, `mid`, `epi` and the `manufactured` sine
fixture. Regional/personalized templates require a separate spatial-gradient
model and are rejected. The existing pseudo-ECG observer retains its regional
path. See the [mathematical and verification guide](../../../../docs/EIKONAL_LEAD_FIELD_ECG.md)
and [maintained dual tutorial](../../../../tutorials/manufacturedSolutions/eikonalLeadFieldECG/README.md).
