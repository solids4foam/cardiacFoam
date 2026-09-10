# listCellModelsVariables

Prints and writes the complete variable inventory of the configured ionic
model and, if present, the active-tension model. Useful for discovering
which state, algebraic, and constant variables are available for
`outputVariables` in `electroProperties`.

## What it does

1. Reads `constant/physicsProperties` to identify the active physics setup.
2. Reads `constant/electroProperties` and constructs the configured ionic model.
3. If an `activeTensionModel` block is present, constructs that model too.
4. Reports three variable groups for each model:
   - **Constants** — fixed model parameters (with initial values)
   - **States** — ODE state variables (with initial values)
   - **Algebraic** — derived quantities computed each step
   - **Rates** — (active-tension models only) rate-of-change variables
5. Writes the same report to `postProcessing/listCellModelsVariables.txt`.

Supported `physicsProperties.type` values:

- `electroModel`
- `electroMechanicalModel`

The utility reads the canonical `myocardiumSolver` selector and constructs the
ionic model from its matching `<myocardiumSolver>Coeffs` subdictionary. The
legacy `electroModel` selector remains supported when it has a matching
`<electroModel>Coeffs` subdictionary.

## Usage

```bash
listCellModelsVariables -case <caseDir>
```

## Output

- Console: full variable report with indexed names and initial values.
- File: `postProcessing/listCellModelsVariables.txt` (same content).

## Notes

- Runs in serial only (`noParallel`).
- The model is instantiated with a single integration point and a nominal
  step size of `0.01` ms; no ODE steps are actually advanced.
- Variable names reported here are the keys accepted by `outputVariables`
  in the solver coefficients dictionary.
