#!/usr/bin/env python3
"""Write the heart conductivity tensors of the case from its fibre field.

Reads the fibre/cross-fibre conductivities of constant/heartConductivity and the
cell field 0/fiber, and writes 0/ConductivityIntracellular and
0/ConductivityExtracellular (volSymmTensorField, S/m):

    G = sigmaCross I + (sigmaFibre - sigmaCross) f f

which equals sigmaFibre f f + sigmaCross (s s + n n) for the orthonormal frame
(fiber, sheet, sheetNormal) of this mesh. Every patch face takes the tensor of
its owner cell (calculated patches). The torso uses bodyAndOrgansConductivity;
the tensors are written on every cell, and only the myocardium reads them.

Run from the case directory (Allrun does it). Standard library only.
"""
import re
import sys


def read(path):
    with open(path) as fh:
        return fh.read()


def header(text):
    return text[: text.index("}", text.index("FoamFile")) + 1]


def scalar_pair(dictText, key):
    m = re.search(rf"^\s*{key}\s+\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)\s*;", dictText, re.M)
    if not m:
        sys.exit(f"heartConductivity.py: constant/heartConductivity has no '{key} (fibre crossFibre);'")
    return float(m.group(1)), float(m.group(2))


def labels(path):
    text = read(path)
    body = text[text.index("}", text.index("FoamFile")) + 1:]
    m = re.search(r"(\d+)\s*\(", body)
    return [int(v) for v in body[m.end():].split(")", 1)[0].split()]


props = read("constant/heartConductivity")
sets = {
    "ConductivityIntracellular": scalar_pair(props, "intracellular"),
    "ConductivityExtracellular": scalar_pair(props, "extracellular"),
}

fiberText = read("0/fiber")
m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", fiberText)
nCells = int(m.group(1))
fib = [tuple(map(float, v)) for v in re.findall(
    r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", fiberText[m.end():])[:nCells]]
if len(fib) != nCells:
    sys.exit("heartConductivity.py: could not read 0/fiber")

owner = labels("constant/polyMesh/owner")
boundary = read("constant/polyMesh/boundary")
patches = [(p.group(1), int(p.group(2)), int(p.group(3))) for p in re.finditer(
    r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", boundary)]


def tensor(f, sf, sc):
    d = sf - sc
    x, y, z = f
    return (sc + d * x * x, d * x * y, d * x * z, sc + d * y * y, d * y * z, sc + d * z * z)


def fmt(t):
    return "(" + " ".join(f"{v:.12g}" for v in t) + ")"


head = header(fiberText)
for field, (sf, sc) in sets.items():
    cells = [tensor(f, sf, sc) for f in fib]
    out = [head.replace("volVectorField", "volSymmTensorField").replace("object      fiber;", f"object      {field};"),
           "\n\ndimensions      [-1 -3 3 0 0 2 0];\n\n",
           f"internalField   nonuniform List<symmTensor>\n{nCells}\n(\n", "\n".join(map(fmt, cells)), "\n)\n;\n\nboundaryField\n{\n"]
    for name, nFaces, start in patches:
        vals = "\n".join(fmt(cells[owner[f]]) for f in range(start, start + nFaces))
        out.append(f"    {name}\n    {{\n        type            calculated;\n"
                   f"        value           nonuniform List<symmTensor>\n{nFaces}\n(\n{vals}\n)\n;\n    }}\n")
    out.append("}\n")
    with open(f"0/{field}", "w") as fh:
        fh.write("".join(out))
    print(f"heartConductivity.py: 0/{field}: fibre {sf} S/m, cross-fibre {sc} S/m, {nCells} cells")
