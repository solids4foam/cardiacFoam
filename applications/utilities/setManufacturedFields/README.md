# setManufacturedFields

Writes the conductivity fields and the cell decomposition that a manufactured-solution verification model prescribes.

```bash
setManufacturedFields -case ./myCase
```

The utility reads `constant/electroProperties`, takes the dictionary of the selected myocardium solver, selects the verification model it names, and calls that model's field writer on the case mesh. The verification model library must be listed under `libs` in `system/controlDict`.

It writes `0/ConductivityIntracellular`, `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity`, and `constant/cellDecomposition` for `method manual` decomposition. Existing files with these names are overwritten.

When an ECG domain of that solver has the verification model `manufacturedLeadFieldECGVerifier` with `reference sphericalShell`, the utility also writes the three conductivity fields of that solution (Gi, Ge and the bath conductivity at the cell centres); the myocardium then needs no verification model of its own.

The utility is serial: run it on the undecomposed case, before `decomposePar`. Running it with `-parallel` is a fatal error.

Failure modes:

- no verification model type, or `none`, and no sphericalShell lead-field ECG domain: the utility stops and names the missing dictionary;
- the verification model library is not loaded: the model table is empty and the error says so;
- the selected model does not write manufactured fields: the model's own error is raised.
