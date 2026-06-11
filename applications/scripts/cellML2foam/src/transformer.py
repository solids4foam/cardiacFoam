#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     transformer
#
# Description
#     Transforms generated C source code for OpenFOAM compatibility.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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
