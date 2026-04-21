# modules/physicsModel architecture

This module provides a lightweight fallback `physicsModel` implementation so
`cardiacFoam` can run electrophysiology workflows without a full `solids4foam`
installation.

## Why this exists

`cardiacFoam` uses `physicsModel::New(...)` as the top-level runtime entrypoint.
In full multi-physics setups this comes from `solids4foam`; in electro-only setups
this module provides the minimum required API.

## Directory structure

```text
modules/physicsModel/
├── src/solids4FoamModels/physicsModel/physicsModel.H/.C
├── src/solids4FoamModels/Make/
├── etc/wmake-options
└── README.md
```

## Runtime behavior

- Reads `constant/physicsProperties`.
- Selects model type from key `type`.
- Constructs runtime-selected model through the physicsModel selection table.
- Provides shared top-level hooks used by `cardiacFoam`:
  - `setDeltaT(...)`
  - `evolve()`
  - `updateTotalFields()`
  - `writeFields(...)`
  - `end()`

`electroModel` registers into this selection system, so electro workflows remain
consistent with the broader solids4foam architecture.

## Scope

This is intentionally minimal and focused on enabling electro runs plus build
portability. Full solid/FSI functionality still requires full `solids4foam`.
