# leadFieldECG

Passive ECG observation for `monodomainSolver` and `bidomainSolver`.

```foam
ecgDomains
{
    ECG { ecgSolver leadFieldECG; }
}
```

Reads `constant/leadFieldProperties` and the precomputed fields created by
`buildLeadFields`. Physical Gi/Ge are mandatory through the existing explicit
`conductivitySource field|uniform` contract. A separate observation
configuration prevents monodomain effective conductivity from being interpreted
as Gi. Weighted fields already contain Gi and reciprocal-current normalization.

`leadFields` maps ordered channel names to native vector fields in
`fieldInstance` (default "0"). Output is `postProcessing/leadFieldECG.dat`,
with time in seconds and channel voltages in volts. Channels are lead
combinations, not fictitious electrode positions. Without `sampling`, output
uses the normal case write times. Optional `sampling { start ...; end ...;
deltaT ...; }` emits actual step times on its independent schedule; deltaT must
be at least the solver timestep. No temporal interpolation is implied.

Restart appends only after validating the existing header, finite monotone
samples and that the trace does not extend beyond the restart time. Starting
again at time zero creates a new trace. Static fields reject moving/changing
meshes and mismatched geometry, topology, Gi/Ge or construction settings.
The per-cell fingerprint uses three 16-bit pieces, safe at ordinary ASCII
field precision. All ranks project local cells; only the master writes.

Heart-cellZone runs import fields through the existing subset/base addressing.
Heart-only runs select `importMode heartOnly` and use explicitly exported
fields on their own myocardial mesh; `originalConductor` is the default and
checks the complete support geometry, actual heart zone and conductor schemes.
The heart-only mode retains the historical exterior-conductor provenance
without loading that exterior mesh. Regenerate if that exterior model changes. Arbitrary external mesh transfers
are not automatic: map and validate those fields separately. Keep native mesh
coordinates and physical snapshots lossless during decomposition (use
`writePrecision 17` even with binary fields, because uniform fields have
textual values); a changed geometry requires new
sensitivities.

The manufactured mono/bidomain tutorials verify the source and vector-field
integral independently. Bounded-conductor configuration cases exercise
construction and Wilson weights. Existing pseudo/torso/eikonal selectors and
scientific defaults are preserved; eikonal lead-field support is deferred.

The observation gradient is selected by `gradSchemes.leadFieldGradVm`.
The private observation field copies the source owner's published wall states
under its own field key. Scheme selection and gradient identity remain
separate, allowing `sealedWallLeastSquares` to reconstruct wall cells without
using sealed-face values. The source must publish its policy; absence is an
error rather than an implicit all-stock assumption. The manufactured cases
select `sealedWallLeastSquares`; the orthogonal column uses `Gauss linear`.
Neither source voltage nor a registered source gradient is updated.
The generator rejects stencils that cannot resolve every active geometric
direction instead of silently zeroing missing components.

For a standalone heart case, first construct fields on the full conductor,
then make a separate working copy and subset its mesh and native observation
fields with `subsetMesh -overwrite -patch <non-empty-patch> -zone myocardium`.
Select `importMode heartOnly` in ECG and remove the primary solver's cellZone
selection when its whole mesh is now myocardium. Keep the observation
properties and metadata as the historical conductor definition. A cut boundary
must retain the intended physical condition; do not turn it into an empty
patch. The configuration-check export hook demonstrates this workflow.

Restart append is verified on this standalone workflow. The PR59 source
reader maps base-mesh restart values into its active cellZone. One ECG domain
of this solver type should own `leadFieldECG.dat`, as with the existing
solver-based output names.
