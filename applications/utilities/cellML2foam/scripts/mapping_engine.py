import os
import re
import sys
from collections import OrderedDict

DEBUG = False

# Extra constants that are not part of updateConstants()
EXTRA_CONSTANTS = [
    "stim_start",
    "stim_duration",
    "stim_amplitude",
    "stim_period_S1",
    "nstim1",
    "stim_period_S2",
    "nstim2",
]
EXTRA_ALGEBRAIC = [
    "Istim",
    "Iion_cm",
]


# ============================================================
# Utilities
# ============================================================

def debug(msg: str, verbose: bool):
    if verbose:
        print(msg)


def parse_model_from_output(output_c):
    """
    Expect output file name: <ModelName>_<Year>.H
    """
    base = os.path.basename(output_c)

    if not base.endswith(".H"):
        raise ValueError("Output file must end with .H")

    stem = base[:-2]  # remove .H

    if "_" not in stem:
        raise ValueError("Output file name must be <ModelName>_<Year>.H")

    model_name, year_str = stem.rsplit("_", 1)

    if not year_str.isdigit():
        raise ValueError("Year must be numeric")

    return model_name, int(year_str)


# ============================================================
# Stimulus helper
# ============================================================

def emit_compute_istim(f):
    f.write(
        """
inline Foam::scalar computeIstim(Foam::scalar t, const double* C)
{
    return Foam::stimulusIO::computeStimulus
    (
        t,
        C[stim_start],
        C[stim_period_S1],
        C[stim_duration],
        C[stim_amplitude],
        Foam::label(C[nstim1]),
        C[stim_period_S2],
        Foam::label(C[nstim2])
    );
}
"""
    )


def emit_dependency_map(f, model_name):
    f.write(
        f"""
static inline const Foam::HashTable<Foam::wordList>&
{model_name}DependencyMap()
{{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {{
        }}
    );

    return dep;
}}
"""
    )


# ============================================================
# Step 1: Load numeric → symbolic mapping
# ============================================================

def load_state_map(fname="state_map.txt"):
    mapping = {}
    with open(fname) as f:
        for line in f:
            if line.strip():
                idx, name = line.split()
                mapping[int(idx)] = name
    return mapping


def states_from_map(mapping):
    return [mapping[i] for i in sorted(mapping.keys())]


# Step 2: Rewrite numeric STATES/RATES
# ============================================================

def rewrite_states_rates(src, mapping):
    pattern = re.compile(r"\b(STATES|RATES)\[(\d+)\]")

    def repl(m):
        arr, idx = m.groups()
        idx = int(idx)
        if idx not in mapping:
            raise RuntimeError(f"Unmapped {arr}[{idx}]")
        return f"{arr}[{mapping[idx]}]"

    out = pattern.sub(repl, src)

    if re.search(r"\b(STATES|RATES)\[\d+\]", out):
        raise RuntimeError("Unreplaced numeric STATES/RATES found")

    return out


# Step 3: Extract symbolic names
# ============================================================

def extract_names(src, array):
    pattern = re.compile(rf"\b{array}\[([A-Za-z_][A-Za-z0-9_]*)\]")
    names = OrderedDict()
    for m in pattern.finditer(src):
        names[m.group(1)] = None
    return list(names.keys())


# Step 4: Extract function blocks
# ============================================================
def extract_function(src, func_name):
    """
    Extract a full C function starting from 'void <func_name>'
    and ending at the matching closing brace.
    """

    # Look for 'void' followed by the function name
    pattern = re.compile(
        rf"void\s*\n?\s*{re.escape(func_name)}\s*\(",
        re.MULTILINE
    )

    m = pattern.search(src)
    if not m:
        raise ValueError(f"Function '{func_name}' not found")

    idx = m.start()  # preserves 'void'

    start = src.find("{", idx)
    if start == -1:
        raise ValueError(f"Opening brace not found for '{func_name}'")

    brace = 0
    for i in range(start, len(src)):
        if src[i] == "{":
            brace += 1
        elif src[i] == "}":
            brace -= 1
            if brace == 0:
                return src[idx:i + 1]

    raise ValueError(f"Unmatched braces in '{func_name}'")


def extract_body_only(func_text):
    return func_text[func_text.find("{") + 1:func_text.rfind("}")].strip()


# Step 5: Extract CONSTANTS from updateConstants()
# ============================================================

def extract_constants_from_update_constants(src):
    block = extract_function(src, "updateConstants")
    pattern = re.compile(r"\bCONSTANTS\[([A-Za-z_][A-Za-z0-9_]*)\]")
    names = OrderedDict()

    for m in pattern.finditer(block):
        names[m.group(1)] = None

    for c in EXTRA_CONSTANTS:
        names.setdefault(c, None)

    return list(names.keys())


# Step 6: Rewrite function signatures
# ============================================================

def rewrite_init_consts_signature(src, model_id):
    pattern = re.compile(
        r"void\s*\n?\s*updateConstants\s*\([^)]*\)\s*\{",
        re.MULTILINE,
    )

    replacement = (
        f"void\n"
        f"{model_id}initConsts("
        f"double* CONSTANTS, double* RATES, double* STATES, "
        f"int tissueFlag, const Foam::dictionary& stimulus"
        f")\n{{"
    )

    return pattern.sub(replacement, src, count=1)


def rewrite_compute_variables_signature(src, model_id):
    pattern = re.compile(
        r"static\s+int\s+rhs\s*\([^)]*\)\s*\{",
        re.MULTILINE,
    )

    replacement = (
        f"void\n"
        f"{model_id}computeVariables("
        f"double VOI, double* CONSTANTS, double* RATES, "
        f"double* STATES, double* ALGEBRAIC, "
        f"int tissueFlag, bool solveVmWithinODESolver"
        f")\n{{"
    )

    return pattern.sub(replacement, src, count=1)


def normalize_void_signature(func_block):
    s = func_block.lstrip()
    if s.startswith("void\n"):
        return func_block
    return re.sub(r"^\s*void\s+([A-Za-z_]\w*)\s*\(",
                  r"void\n\1(", func_block)


# Step 7: Header helpers
# ============================================================

def emit_enum(f, name, symbols):
    base = name.replace("_INDEX", "")
    f.write(f"enum {name} {{\n")
    for s in symbols:
        f.write(f"    {s},\n")
    f.write(f"    NUM_{base}\n}};\n\n")


def emit_names_array(f, array, enum, symbols):
    f.write(f"static const char* {array}[NUM_{enum}] = {{\n")
    for s in symbols:
        f.write(f'    "{s}",\n')
    f.write("};\n\n")


def emit_function_prototypes(f, model_id):
    f.write(
        f"""
void
{model_id}initConsts(double* CONSTANTS,double* RATES,double* STATES,int tissueFlag,const Foam::dictionary& stimulus);

void
{model_id}computeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC,int tissueFlag,bool solveVmWithinODESolver);
"""
    )


def emit_openfoam_algebraic_tail():
    return """

    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[0] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }
"""


def generate_header(fname, states, algebraic, constants, model_name):
    with open(fname, "w") as f:
        f.write("#pragma once\n\n")
        emit_enum(f, "STATES_INDEX", states)
        emit_enum(f, "ALGEBRAIC_INDEX", algebraic)
        emit_enum(f, "CONSTANTS_INDEX", constants)
        emit_function_prototypes(f, model_name)


# Final C assembly
# ============================================================
def build_final_c(init_block, computeVariables_block, init_values_block):

    # Ensure correct headers
    init_block = normalize_void_signature(init_block)
    computeVariables_block = normalize_void_signature(computeVariables_block)

    init_body = extract_body_only(init_block)
    init_vals = extract_body_only(init_values_block)

    # Reconstruct init function explicitly
    init_header = init_block.split("{", 1)[0].rstrip() + "\n{"

    final_init = (
        init_header
        + "\n"
        + init_body
        + "\n\n    /* --- Initial values --- */\n\n"
        + init_vals
        + "\n}\n"
    )

    compute_header = computeVariables_block.split("{", 1)[0].rstrip() + "\n{"
    compute_body = extract_body_only(computeVariables_block)
    # Remove Myokit rhs return
    compute_body = re.sub(r"\breturn\s+0\s*;", "", compute_body)

    final_compute = (
        compute_header
        + "\n"
        + compute_body
        + emit_openfoam_algebraic_tail()
        + "\n}\n"
    )

    return final_init + "\n\n" + final_compute


def run_mapping(input_c, output_c, verbose: bool = False):
    global DEBUG
    DEBUG = verbose
    model_name, year = parse_model_from_output(output_c)
    model_id = f"{model_name}_{year}"
    output_h = f"{model_id}Names.H"

    with open(input_c) as f:
        src = f.read()

    mapping = load_state_map()

    rewritten = rewrite_states_rates(src, mapping)
    rewritten = rewrite_init_consts_signature(rewritten, model_name)
    rewritten = rewrite_compute_variables_signature(rewritten, model_name)

    states = states_from_map(mapping)
    algebraic = extract_names(rewritten, "ALGEBRAIC")
    for name in EXTRA_ALGEBRAIC:
        if name not in algebraic:
            algebraic.append(name)
    constants = extract_constants_from_update_constants(src)

    init_block = extract_function(rewritten, f"{model_name}initConsts")
    computeVariables_block = extract_function(
        rewritten, f"{model_name}computeVariables"
    )
    init_vals = extract_function(rewritten, "default_initial_values")

    debug(computeVariables_block, verbose)

    final_c = build_final_c(init_block, computeVariables_block, init_vals)

    with open(output_c, "w") as f:
        f.write(f'#include "{output_h}"\n')
        f.write('#include "stimulusIO.H"\n\n')
        emit_names_array(f, f"{model_name}STATES_NAMES", "STATES", states)
        emit_names_array(f, f"{model_name}ALGEBRAIC_NAMES", "ALGEBRAIC", algebraic)
        emit_dependency_map(f, model_name)
        emit_compute_istim(f)
        f.write(final_c)

    generate_header(output_h, states, algebraic, constants, model_name)

    debug(f"Wrote: {output_c}", verbose)
    debug(f"Wrote: {output_h}", verbose)

    return {
        "model": model_name,
        "year": year,
        "header": output_h,
        "source": output_c,
    }


def main():
    input_c, output_c = sys.argv[1:]
    run_mapping(input_c, output_c, verbose=True)


if __name__ == "__main__":
    main()
