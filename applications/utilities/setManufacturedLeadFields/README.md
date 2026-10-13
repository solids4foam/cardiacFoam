# setManufacturedLeadFields

Verification fixture generator for the unit-cube FDA voltage family. Reads
explicit physical conductivities and `leadFields` from
`constant/leadFieldProperties`; `manufacturedScales` supplies one signed scale
per channel. Writes conductivity-weighted vectors with frequencies (1,2,3)
and the same mesh/conductivity provenance used by the production importer.

Run after `blockMesh`, before `decomposePar` or `cardiacFoam`. The fields are
synthetic integration fixtures, not clinical ECG sensitivity fields. The
manufactured mono/bidomain lead-field tutorials supply the complete workflow.

`-checkSealedWall` generates constant sensitivity fixtures and verifies affine
voltage reconstruction on sealed walls, with a deliberately stale registered
source gradient. It publishes the source's wall states and exercises the
private observation field's own `sealedWallLeastSquares` identity. Optional
extracellular source snapshots are supplied, but are not used to close an
observation offset: both source Vm and phiE arrays must remain unchanged.
Allrun restores the cosine-product fixtures after this focused check.

`-eikonal` writes constant sensitivities `(1,2,3)` and its signed/scaled
variants, then exercises the production eikonal lead-field observer with
four uniform templates. It checks an affine activation field, compiled-table
units/slopes, an activation-only provider, a stale source gradient, unchanged
source arrays and repeated solves. The
[eikonal dual tutorial](../../../tutorials/manufacturedSolutions/eikonalLeadFieldECG/README.md)
uses this option before the PDE run.

`-eikonalGuard regionalTemplate|unknownTemplate|invalidSampling` exercises
constructor rejection without changing fields. This option must fail with
the corresponding diagnostic and is used by the maintained regression.
