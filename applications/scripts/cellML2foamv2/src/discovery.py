"""
discovery.py - Automated discovery of cardiac model components from Myokit.
"""

def find_variable_by_label(model, label):
    """Finds a variable in a Myokit model by its semantic label."""
    try:
        return model.label(label)
    except KeyError:
        return None

def find_voltage(model):
    """Heuristic logic to find the membrane voltage variable."""
    # 1. Check label (Strongest hint)
    v = find_variable_by_label(model, 'membrane_voltage')
    if v:
        return v
    
    # 2. Check common names in state variables
    for s in model.states():
        name = s.name().lower()
        if name in ['v', 'vm', 'v_m', 'membrane_v', 'membrane_voltage']:
            return s
            
    # 3. Fallback: first state variable (usually V in cardiac models)
    states = list(model.states())
    if states:
        return states[0]
        
    return None

def find_stimulus(model):
    """Heuristic logic to find the stimulus current variable."""
    # 1. Check labels
    labels = ['membrane_stimulus_current', 'stimulus_protocol', 'pace']
    for label in labels:
        v = find_variable_by_label(model, label)
        if v:
            return v
            
    # 2. Check common names
    for v in model.variables():
        name = v.name().lower()
        if 'istim' in name or 'i_stim' in name or 'pace' == name:
            return v
            
    return None

def find_ionic_current(model, voltage_var):
    """
    Analyzes the derivative of voltage to find the ionic current component.
    """
    if not voltage_var:
        return None
        
    # 1. Check for explicit label
    i_ion = find_variable_by_label(model, 'membrane_ionic_current')
    if i_ion:
        return i_ion
        
    # 2. Structural heuristics on dV/dt
    # dot(V) = -(Iion + Istim) / Cm
    try:
        rhs = voltage_var.rhs()
        # Find all variables involved in the RHS of dot(V)
        refs = list(rhs.references())
        
        # Candidate ionic current is often named i_ion, i_total, i_all, etc.
        for ref in refs:
            name = ref.name().lower()
            if any(term in name for term in ['i_ion', 'i_total', 'i_all', 'i_total_ion']):
                return ref
    except:
        pass
        
    return None

def discover_all(model):
    """
    Performs discovery on the model and returns a mapping of 
    OpenFOAM standard keys to Myokit variable names.
    """
    results = {
        'Vm': None,
        'Istim': None,
        'Iion': None,
    }
    
    if not model:
        return results
        
    v = find_voltage(model)
    if v:
        results['Vm'] = v.name().replace('.', '_')
        
    s = find_stimulus(model)
    if s:
        results['Istim'] = s.name().replace('.', '_')
        
    i = find_ionic_current(model, v)
    if i:
        results['Iion'] = i.name().replace('.', '_')
        
    return results
