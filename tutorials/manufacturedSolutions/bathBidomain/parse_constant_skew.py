import sys
import re

path = "constant/skewness"
try:
    with open(path) as f: text = f.read()
except:
    path = "constant/polyMesh/skewness"
    with open(path) as f: text = f.read()

m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", text)
if m:
    vals = [float(m.group(1))]
else:
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\(", text)
    if m:
        n = int(m.group(1))
        start = text.index("(", m.end() - 1) + 1
        end = text.index(")", start)
        vals = [float(v) for v in text[start:end].split()]

print(f"Mean: {sum(vals)/len(vals):.6f}, Max: {max(vals):.6f}")
