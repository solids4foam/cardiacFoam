# leadFieldECG

The lead-field evaluation of `torsoECG`: the ECG of a heart+torso passive
conductor from lead fields precomputed by `buildLeadFields`. Each channel is
the sum, over the heart cells, of the electrode's lead field times the heart
source of the current `Vm`, plus the constant term of the conductor's
boundary data and reference. No elliptic problem is solved per step, and
one build serves every run on the same heart, torso, conductivities and
electrodes. The boundary data are therefore constant: a surface current
that varies in time stops the run (`passiveConductor`); the bath
(`torsoCoupling bath`) takes it.

- `leadFieldConductor` builds the heart+torso conductor (`passiveConductor`,
  shared with the bath) for an ECG domain, with its electrodes and the
  signature that ties the lead fields to the case.
- `leadFieldAdjoint` transposes the explicit part of the conductor's
  operator for `buildLeadFields`.
- `leadFieldECG` loads and evaluates the lead fields. At load the signature
  must match the case, and each lead field's fingerprint (sums of `z V f`
  over the cells) must match the one `buildLeadFields` recorded in
  `leadFieldBuildRecord`. A renumbered or redistributed mesh keeps the
  signature's totals but fails the fingerprint.

The lead field of electrode `e` is the discrete adjoint `z_e = A^-T e_e` of
the finite-volume operator `A` of the conductor, with its explicit
non-orthogonal, cross-diffusion and interface corrections transposed. The
channel then equals the electrode value of the converged forward solve
`A phi = b` of the same conductor on any mesh; `A` is not symmetric on
general meshes, so the lead field of `A` itself would not.

Electrodes lie in the torso: an electrode outside the mesh or in a
myocardium cell stops the run with a fatal error. A bidomain myocardium
whose heart potential takes the bath's trace on the heart wall is not
supported: the conductor's heart companion is zeroGradient.

## The heart wall

The forward problem is the extracellular problem for a given `Vm`:
`div((Gi + Ge) grad phi) = -div(Gi grad Vm)` in the heart,
`div(sigma grad phi) = 0` in the torso, and `phi` and the total current
continuous at the heart wall. The heart side of the heart-torso faces
carries `Gi + Ge`, and the heart source moves the wall current
`q = n.Gi grad(Vm)` of each wall face to the torso cell
(`passiveConductor::transferWallLayer`).

The bath (`torsoCoupling bath`) solves the same problem but writes the wall
with `Ge`. It solves `Vm` and `phi` together, so the pair satisfies
intracellular insulation and only extracellular current crosses the wall.
A `Vm` given to the lead fields does not satisfy that identity with the
torso's `phi`, so only the total-current form holds. For the bath's own `Vm`
the two forms are the same equation. See `passiveConductor.H`.
