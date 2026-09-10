# setTorsoOrganConductivityField

Create a static torso/bath conductivity field from named mesh `cellZones`.

Default dictionary:

```text
system/setTorsoOrganConductivityFieldDict
```

Example:

```foam
fieldName bodyAndOrgansConductivity;

cellZones
{
    torso  0.2;
    lungs  0.04;
    blood  0.7;
    bone   0.02;
}

// Optional. Only used for cells not covered by the zones above.
defaultSigma 0.2;
```

The utility fails if a configured cellZone name is not present in the mesh. It
also fails if cells remain uncovered and no `defaultSigma` is supplied.

`bodyAndOrgansConductivity` is also the default field read by
`bidomainSolverCoeffs.bathPotentialDomain`. Use
`bathConductivityField <name>;` inside that block if a case needs a different
field name.
