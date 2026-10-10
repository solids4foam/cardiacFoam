# Transmission manufactured-solution generator

`transmission_mms.py` builds, with sympy, the exact solution of the
bath-bidomain transmission problem on the oblique slab: a heart slab between
two bath slabs with planar interfaces oblique to the mesh axes, rotated
anisotropic heart conductivities, spatially varying bath conductivities, and
fields that satisfy potential continuity, extracellular flux continuity and
intracellular sealing on both interfaces. It checks the construction (interface
conditions, sources against finite differences, boundedness, periodicity,
richness of the tangential conormal flux, common-mode invariance of Vm, the
affine solution, the ionic closure) and writes the C++ reference header
`../manufacturedTransmissionBathBidomainReference.H` used by the transmission
verifier.

Each entry of `selfTestSamples` holds a function index, a parameter set
(0 `defaultParams()`, 1 `alternateParams()`), the solution (0 affine,
1 levelSet), the zone, a component and the reference value. For
`interfaceFlux` the component `c` is evaluated with the unit normal along
axis `c`. Functions without a time argument (`Gi`, `Ge`, `sigmaB`, `G`) are
sampled at the first time only; the others at two times.

`--write` runs every generator check first and writes nothing if one fails.
A relative `--write` or `--check` path is taken against this directory.

Regenerate the header (from this directory):

```bash
python transmission_mms.py --write ../manufacturedTransmissionBathBidomainReference.H
```

Check that the tracked header matches the generator (exit status 1 on a
difference):

```bash
python transmission_mms.py --check ../manufacturedTransmissionBathBidomainReference.H
```

Run the generator checks:

```bash
python test_transmission_mms.py
```

Requires sympy, mpmath and numpy; the header records the sympy version that
produced it. `--check` and the header comparison in `test_transmission_mms.py`
compare text exactly, including that version line and the common
subexpression output, so they hold only for the recorded sympy version.
