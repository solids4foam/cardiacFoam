# bathBidomainInterfaceMetrics

Compute manufactured heart, bath, interface, and exterior-boundary diagnostics
from a serial or reconstructed bath-bidomain result.

```bash
bathBidomainInterfaceMetrics -latestTime
```

The case must contain `phiE`, `VmGlobal`, and `sigmaTotal` at the selected time plus
`myocardium` and `bath` cellZones. Output is written to:

```text
postProcessing/bathBidomainInterfaceMetrics.csv
```

Cell norms are volume weighted and interface/boundary norms are area weighted.
One-sided interface fluxes are reconstructed independently from cell gradients;
the reported flux jump therefore tests constitutive reconstruction rather than
merely restating the conservative assembled face flux.
The utility also reconstructs the intracellular interface leakage
`Gi*grad(Vm + phiE).n`, whose manufactured value is zero.

The assembled current density is the face flux of
`fvm::laplacian(Sef, phiE)`, with `Sef` the face conductivity of the
heart+bath conductor (`passiveConductor`) and the standard cell gradient of
`phiE` in its corrections: the flux of the `distanceWeightedHarmonic`
operator. Under `conormalHarmonic` the solver's corrections use the
jump-aware gradient, which the utility does not apply, so the CSV carries no
interface option.
