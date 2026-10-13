# buildLeadFields

Constructs static ECG sensitivity vectors on a fixed passive conductor mesh.
Run after mesh/conductivity preparation, before `decomposePar` or `cardiacFoam`.
The standard OpenFOAM `-parallel` option also supports construction on an
already decomposed conductor. No cardiac voltage is advanced or modified.

Configuration is `constant/leadFieldProperties`. `conductivitySource` must be
`uniform` (explicit dimensioned `conductivityIntracellular` and
`conductivityExtracellular`) or `field` (canonical `ConductivityIntracellular`
and `ConductivityExtracellular`, with the existing legacy names accepted).
These are physical conductivities, independent of monodomain effective
propagation conductivity. They use S/m and native symmetric tensors.

The current generator requires myocardium fully embedded in a passive
conductor; exposed myocardium on a physical exterior patch is rejected.
Lower-dimensional empty patches and coupled patches retain their constraint
semantics. Supply `cellZone myocardium` and a scalar
`bodyAndOrgansConductivity` field (or `bathConductivityField` override) on the
base mesh, covering all passive cells with positive S/m values. Within the
heart the operator combines Gi and Ge. Heart/conductor interface faces use
only extracellular conductivity with distance-weighted harmonic interpolation;
intracellular current does not cross those faces. Physical outer patches use
ordinary zeroGradient insulation for the scalar passive exterior; coupled/empty patch types are retained.

Supply `electrodePositions`, balanced `leads` dictionaries and a matching
`leadFields` mapping. `phiERefPoint` must lie in the mesh. Electrode locations
must also lie in the mesh; outside/skin locations can use an explicit positive
`maximumElectrodeDistance` in metres for a bounded nearest-cell approximation.
Reciprocal excitation and recording use the same owning-cell observation.

Example of Wilson's reference (all electrode coordinates must be supplied):

```foam
leads
{
    V1 { E1 1; RA -0.333333333333333333; LA -0.333333333333333333; LL -0.333333333333333333; }
}
leadFields { V1 leadField_V1; }
fieldInstance "0";
```

Use the same construction for E2-E6/V2-V6. The electrode weights set the lead's
polarity and scale. The reciprocal problem is normalized per ampere and the
output is `L=-Gi^T grad(z)` in m^-2, excluding cell volume. Production projection
is `sum(volume * L dot grad(Vm))`, yielding volts. Do not multiply by Gi or the
pseudo-ECG surrounding-medium normalization again.

Configure a `leadPotential` solver in `system/fvSolution`, and the normal
Laplacian/gradient schemes in `fvSchemes`. `maxCorrectors` (default 100) and
`correctionTolerance` (default 1e-10) control converged corrected assemblies.
A failed correction convergence stops construction. The gauge is zero at
`phiERefPoint`; balanced lead differences are gauge independent.

Myocardial gradients are reconstructed with least squares from myocardial
neighbours, including processor neighbours. This preserves nonzero reciprocal
current at the heart interface and avoids mixing the bath-side slope into the
heart. No artificial no-flux condition is imposed at the restriction boundary.

Outputs include the vectors, optional reciprocal potentials (`writePotentials
yes`, with their required auxiliary fields), physical conductivity
snapshots, geometry/topology fingerprints and `constant/leadFieldMetadata`.
The runtime importer checks these before projection, including original
construction settings. Preserve lossless mesh coordinates when transferring or
decomposing anatomical meshes. The bounded column configuration checks provide
a complete runnable example for both monodomain and bidomain.

This implementation uses a finite-volume reciprocal field and cell-gradient
quadrature. It is not the full discrete transpose of cardiacFoam's corrected
coupled bath assembly; exact finite-grid agreement with that solve is not
claimed. Validate spatial convergence and the cardiac source model separately.

[examples/wilsonLeads](examples/wilsonLeads) supplies all six E1-E6/V1-V6 combinations.

The myocardial least-squares normal matrix is checked with SVD (relative rank
cutoff 1e-12). Full-rank gradients use the ordinary inverse; rank-one/two
stencils use the SVD pseudoinverse over their resolved range. This avoids
OpenFOAM cubic eigenvalue/eigenvector failures for rotated repeated eigenpairs.
`buildLeadFields -checkReconstruction` runs four small affine algebra controls
(isotropic, rotated repeated full rank, and rotated rank two/one) before any
conductivity solve. The parallel construction hook requires these checks.
Provenance labels this construction `myocardialNeighbourLeastSquaresSVD`.
