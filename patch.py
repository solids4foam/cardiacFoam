def _single_cell_trace_variables(ionic_model: str | None) -> tuple[str, ...]:
    if ionic_model is None:
        return ()
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    entry = IONIC_MODEL_CATALOG.get(ionic_model)
    if entry is None:
        return ()
    
    # The C++ trace is [Time STATES ALGEBRAIC RATES]
    # We output all state variables, their rates, and algebraics.
    vars_list = ["Time"]
    for s in entry.states:
        vars_list.append(s)
        vars_list.append(f"{s}Rate")
    vars_list.extend(entry.algebraic)
    return tuple(vars_list)
