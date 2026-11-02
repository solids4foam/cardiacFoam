# setManufacturedFields

Writes the conductivity fields and the cell decomposition that a manufactured-solution verification model prescribes.

```bash
setManufacturedFields -case ./myCase
```

The utility reads `constant/electroProperties`, takes the dictionary of the selected myocardium solver, selects the verification model it names, and calls that model's field writer on the case mesh. When the myocardium has no verification model (no `type`, or `none`), it selects the verification model of the only `ecgDomains` entry that has one and calls that model's field writer instead (for example `manufacturedLeadFieldECGVerifier` with `reference sphericalShell`, which writes Gi, Ge and the bath conductivity of its solution at the cell centres). The verification model library must be listed under `libs` in `system/controlDict`.

It writes `0/ConductivityIntracellular`, `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity`, and `constant/cellDecomposition` for `method manual` decomposition. Existing files with these names are overwritten. One model writes the fields: with a myocardium verification model, the ECG domains' models write none.

The utility is serial: run it on the undecomposed case, before `decomposePar`. Running it with `-parallel` is a fatal error.

Failure modes:

- no verification model type, or `none`, in the myocardium and in every ECG domain: the utility stops and names the missing dictionary;
- no myocardium verification model and more than one ECG domain with one: the utility stops and names the domains;
- the verification model library is not loaded: the model table is empty and the error says so;
- the selected model does not write manufactured fields: the model's own error is raised.
