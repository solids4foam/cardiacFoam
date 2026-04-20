import re

def transform_sim_c(src: str) -> str:
    """
    Replaces Coccinelle transformations with pure Python regex.
    Goal: Transform Myokit ANSI C export to OpenFOAM-style C++ logic.
    """

    # 1. STATES: NV_Ith_S(y, i) -> STATES[i]
    src = re.sub(r"NV_Ith_S\s*\(\s*y\s*,\s*(\d+|\w+)\s*\)", r"STATES[\1]", src)

    # 2. RATES: NV_Ith_S(ydot, i) -> RATES[i]
    src = re.sub(r"NV_Ith_S\s*\(\s*ydot\s*,\s*(\d+|\w+)\s*\)", r"RATES[\1]", src)

    # 3. CONSTANTS: AC_... -> CONSTANTS[AC_...]
    # We must be careful not to match #defines or declarations
    # Look for AC_ followed by name, not preceded by #define or type
    src = re.sub(r"(?<!#define\s)(?<!const\s)(?<!double\s)\b(AC_[A-Za-z0-9_]+)\b", r"CONSTANTS[\1]", src)

    # 4. ALGEBRAIC: AV_... -> ALGEBRAIC[AV_...]
    src = re.sub(r"(?<!#define\s)(?<!const\s)(?<!double\s)\b(AV_[A-Za-z0-9_]+)\b", r"ALGEBRAIC[\1]", src)

    # 5. Correct index of STATES/RATES if needed (Myokit often uses y[0], y[1] as standard aliases)
    # But in ANSI C export, it usually maintains NV_Ith_S unless flags are set.

    return src


def run_transformation(file_path, verbose=False):
    if verbose:
        print(f"    Transforming {file_path} (Python-driven rewrite)...")

    with open(file_path, "r") as f:
        src = f.read()

    transformed = transform_sim_c(src)

    with open(file_path, "w") as f:
        f.write(transformed)
