# setManufacturedFields

Writes the conductivity fields and the cell decomposition that a manufactured-solution verification model prescribes.

```bash
setManufacturedFields -case ./myCase
```

The utility reads `constant/electroProperties`, takes the dictionary of the selected myocardium solver, selects the verification model it names, and calls that model's field writer on the case mesh. The verification model library must be listed under `libs` in `system/controlDict`.

It writes `0/ConductivityIntracellular`, `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity`, and `constant/cellDecomposition` for `method manual` decomposition. Existing files with these names are overwritten.

The utility is serial: run it on the undecomposed case, before `decomposePar`. Running it with `-parallel` is a fatal error.

Failure modes:

- no verification model type, or `none`: the utility stops and names the missing dictionary;
- the verification model library is not loaded: the model table is empty and the error says so;
- the selected model does not write manufactured fields: the model's own error is raised.
