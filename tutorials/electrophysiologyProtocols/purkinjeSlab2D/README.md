# purkinjeSlab2D

A small 2-D tissue slab (thin `frontAndBack`-empty block, matching the
`rotorInstability` mesh style) coupled to a minimal 2-node Purkinje network
(`restitutionEikonalSolver1D`, root node 0, one PVJ terminal at node 1) via
bidirectional PVJ coupling.

## Purpose

This is a general-purpose testbed for exercising myocardium/network/coupler
combinations cheaply (small mesh, short runs), not a single fixed scenario.
It currently supports:

- `solver=monodomain` (default) — `monodomainSolver` + BuenoOrovio tissue,
  `eikonalMonodomainPvjCoupler`, `couplingMode bidirectional`. Runs as a real
  multi-step timeline. This is the variant used to investigate retrograde
  tissue -> network activation (the `restitutionEikonalSolver1D` review fix).
- `solver=eikonal` — `eikonalSolver` tissue, `eikonalPvjCoupler`,
  `couplingMode bidirectional`. Kept for comparison only:
  `eikonalMyocardiumDomain::applyModelTimeControls` collapses the whole run
  to a single dimensionless step, so this variant cannot produce a
  multi-step retrograde timeline (confirmed empirically, not just from
  reading the code).

Run with:

```bash
./Allrun                  # solver=monodomain
./Allrun solver=eikonal
```

## Status / open findings

This tutorial surfaced two real issues in the coupling code, found by running
it, not by reading the code:

1. **Fixed**: `myocardiumDomain::activationTime_` (monodomain tissue)
   defaulted to `dimensionedScalar(0.0)` instead of the `-1` "unactivated"
   sentinel used everywhere else in this codebase (`eikonalMyocardiumDomain`,
   `conductionSystemDomain::setTerminalActivationTime`,
   `pvjMapper::gatherActivationTimes`). Under bidirectional PVJ coupling this
   made every not-yet-activated cell look like a genuine activation at
   `t=0`, spuriously firing the network on the very first coupling exchange.
   Fixed in `myocardiumDomain.C` (default changed to `-1`); verified by
   re-running this tutorial with the fix and confirming the PVJ node stays
   at `-1` instead of jumping to `0.0`.

2. **Open**: with the fix applied, the PVJ node correctly stays at `-1` for
   the whole run — but `activationTime` never records a genuine crossing
   *anywhere* in the mesh, even though `Vm` shows real spatial structure
   consistent with a spreading wave. Not yet root-caused: candidates include
   `activationThreshold_`/tissue resting-value mismatch, or
   `updateActivationTime()`'s old/new `Vm` comparison not seeing a real
   transition under the explicit diffusion solve path. This blocks observing
   a genuine future -> due -> stale retrograde transition end-to-end.

## Files

- `constant/electroProperties.monodomain`, `constant/electroProperties.eikonal`
  — solver variants, copied to `constant/electroProperties` by `Allrun`.
- `constant/purkinjeGraph` — 2-node graph (root=0, PVJ terminal=1).
- `constant/physicsProperties`, `system/*` — standard case setup.
