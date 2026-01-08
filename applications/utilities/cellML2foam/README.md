# cellML2foam

`cellML2foam` converts CellML models into OpenFOAM-ready ionic model code. It
is designed for the cardiacFoam workflow and expects the user to create a
`state_map.txt` mapping after inspecting the generated `.mmt`/`ansic` output.

## What it does

- **cellml → mmt** using Myokit
- **mmt → ansic** using Myokit
- **ansic → openfoam** using a state map and Coccinelle rewrites
- Generates model wrapper files from templates

## Requirements

Install the following tools:

- `python3`
- `myokit` (CLI must be on PATH)
- `spatch` (Coccinelle, for source-to-source rewrites)

If `spatch` is unavailable, the pipeline will fail at the ansic → openfoam step.

## How to use

Run the utility from this directory:

```
cd /Users/simaocastro/cardiacFoam/applications/utilities/cellML2foam
```

### 1) Convert CellML to MMT

```
cellML2foam --from cellml --to mmt path/to/model.cellml
```

This writes `model.mmt`.

### 2) Create state_map.txt

Open the `.mmt` or `ansic/sim.c` and list the initial states in order. Create
`state_map.txt` in this directory, for example:

```
0   V
1   CaMKt
2   Nai
3   Nass
```

### 3) Generate OpenFOAM code

```
cellML2foam --from mmt --to openfoam --model ModelName_Year model.mmt
```

Example:

```
cellML2foam --from mmt --to openfoam --model NashPanfilov_2004 nash_panfilov_2004.mmt
```

This produces:

- `<Model>_<Year>.H`
- `<Model>_<Year>Names.H`
- A new model folder with templated wrapper files

## Notes

- The utility expects `state_map.txt` to be in the **current working directory**.
- The `--model` value must be `ModelName_Year` and the year must be numeric.

