# Batched manufactured-monodomain parity

This study uses the existing `monodomainPseudoECG` analytical manufactured
solution to validate `monodomainFDAManufacturedBatched`.

The four JSON specifications copy the scalar Cartesian, temporal,
tetrahedral, and tet timestep-control ladders. They select the batched ionic
model with 25 Euler substeps. Run each specification once with the host-batched
path and once with CUDA visible; generated results stay in separate local
result directories.

`setup/post_processing_manufactured.py` remains the normal convergence
postprocessor: each completed batched case is checked against the same exact
solution and field verifier as its scalar counterpart. The scalar and batched
finite-resolution errors need not be equal, but both must converge to the
same analytical solution. The temporal study establishes whether the batched
ODE time error is sufficiently small for interpreting a spatial slope.

After matching host-batched and CUDA-batched runs complete, use the study
postprocessor for direct backend parity:

```bash
python3 tutorials/manufacturedSolutions/monodomainPseudoECG/setup/studies/monodomainBatchedParity/postprocess_backend_parity.py \
    /path/to/host-batched-case /path/to/cuda-batched-case \
    --times 0.2
```

It compares `Vm`, manufactured fields, ionic current, activation time, and
the saved batched ODE state. The two inputs must use identical mesh, timestep,
coupling scheme, integrator, and substep count.
