# Heart–bath affine interface patch test

This isolated test compiles the current `exposedPhiETrace.C` and calls the
current `extracellularFaceConductivity.H` helper. It does not need an ionic
model or a time step.

The prescribed extracellular potential is continuous and piecewise affine:

- heart, `0 < x < 1`: `phiE = g*y + 0.25*x`;
- left bath, `x < 0`: `phiE = g*y + s_B*x`;
- right bath, `x > 1`: `phiE = g*y + 0.25 + s_B*(x - 1)`;
- `s_B = (0.25*G_e,xx + g*G_e,xy)/sigmaBath`.

The main test uses `g=1`. A `g=0` run is a nontrivial normal-flux
control: the normal gradient still jumps, but the harmonic normal-face
conductivity should reproduce its flux exactly.

Its exact conormal flux is continuous at both interfaces and equals
`0.25*G_e,xx + g*G_e,xy`, while its normal derivative jumps. Consequently, the exact
interior divergence is zero. The test measures the production exposed-face
trace against the exact face value, the flux implied by the production
interpolated tensor against the exact flux, and the OpenFOAM Laplacian
residual in interface-adjacent cells. It uses the conormal bath MMS tensor
and `blockMesh`/`topoSet` files. Samples near the periodic `y` seam are
excluded because this affine probe is deliberately not periodic.

Run from any directory after setting up OpenFOAM:

```bash
wmake applications/utilities/bathInterfacePatchTest
applications/utilities/bathInterfacePatchTest/run.sh 10 20 40
```

`run.sh` fails if the trace, interface flux or adjacent-cell Laplacian
residual exceeds `1e-10` for this exactly affine field. Use
`DIAGNOSTIC_ONLY=1` to print the refinement table without treating a
current inconsistency as a shell failure.

The printed values are diagnostics, not a full cardiac MMS convergence
result. In particular, this aligned Cartesian case does not validate
nonorthogonal or curved interface reconstruction. A physically consistent
affine interface flux should have a flux error near solver roundoff. A
nonzero flux error that stays constant as `N` increases demonstrates an
interface flux inconsistency, irrespective of the observed `Vm` order.
