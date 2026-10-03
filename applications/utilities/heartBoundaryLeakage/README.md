# heartBoundaryLeakage

Net current through the sealed heart boundary.

```bash
heartBoundaryLeakage [-time ...] [-zone myocardium|none] [-exposedPatch <patch>] \
    [-conductivity ConductivityIntracellular] [-fields '(Vm phiE)'] \
    [-refPoint '(x y z)'] [-refValue 0] [-sigmaBath bodyAndOrgansConductivity] \
    [-analyticField xy] [-insulated] [-insulationDelta] \
    [-gradientBias] [-phiETrace zeroGradient|global|conductivityWeighted] \
    [-conductivityExtracellular ConductivityExtracellular] \
    [-exactPhiE bathFDA -mmsSe <s_e> [-mmsK 0.7071067811865476] [-mmsAlpha 0.01]] \
    [-conormalIterations N] [-restDrift] [-restValue -0.084]
```

One line per time, columns in order:

| column | option |
|---|---|
| `time` | |
| `L_<field> Lnorm_<field>` | per field |
| `L` | |
| `delta_<field>=<max>/<scale>@cell<i>(<x>)` | `-insulationDelta`, per field |
| `driftWall driftInner` | `-restDrift` |
| `gradBiasWall gradBiasWallMax gradBiasInner` | `-gradientBias` |
| `globalBiasWall` | `-gradientBias -exactPhiE` |
| `d_r phiRefPredicted phiRefSolved relErr` | `-refPoint` |

`L_<field>` is the sum over heart cells of `V*laplacian(G, field)` with the field `zeroGradient` on every heart patch. `Lnorm_<field>` is the sum of `V*|laplacian(G, field)|`.

- `-insulated` evaluates the Laplacian with `insulatedFaceConductivity` (zero face conductivity on sealed patches).
- `-conormalIterations N` sets `conormalZeroFlux` on every heart patch and runs N patch/gradient passes before the Laplacian.
- `-gradientBias` compares the heart-submesh `grad(phiE)` with a reference: the global-mesh gradient, or with `-exactPhiE bathFDA` the exact gradient of the manufactured bath solution. `-phiETrace` selects the heart `phiE` trace on exposed faces. `conductivityWeighted` weighs the heart and bath cells by `n.G_e.n/d` and `n.sigma.n/d`.
