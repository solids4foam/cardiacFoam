# Electromechanics protocols

Small, fast electromechanics cases, each checking one mechanical ingredient
of the coupled model, such as boundary supports or the active-tension
response.

```text
electromechanicsProtocols/
└── springSupportedSlab/    Niederer slab on spring (solidRobin) end supports
```

All cases here need cardiacFoam built with solids4foam.

## `springSupportedSlab/`

The Niederer et al. (2011) slab, activated by a plane wave along the fibres,
with both fibre-wise ends on `solidRobin` spring supports. Sweeping the
spring stiffness takes the twitch from isometric (no shortening) to free
shortening. The case checks that the spring force computed from the stress
field matches the spring law at every time step.

## Regression

`springSupportedSlab` is wired into `tutorials/Alltest-regression`, as an
expected skip in `lightweight` build mode. See `../README.md` for the full
canonical table.
