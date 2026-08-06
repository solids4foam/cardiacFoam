#!/usr/bin/env python3
"""
make_lbbb_graph.py
==================
Produces a modified purkinjeGraph that simulates Left Bundle Branch Block (LBBB)
by zeroing the conductance of the single bridge edge connecting the His bundle
to the Left Bundle Branch subtree.

Usage
-----
    python3 make_lbbb_graph.py [--dry-run]

Output
------
    constant/purkinjeGraph.lbbb   -- drop-in replacement for constant/purkinjeGraph

Anatomy recap
-------------
  Root node   : 14083  (stimulus site / His bundle entry)
  Bifurcation : 14082  (His bundle → LBB + RBB)
  LBB bridge  : edge 44002 = 4(14082 44015 ...)   <-- zeroed here
  RBB bridge  : edge 43954 = 4(14082 43979 ...)   <-- untouched

Identification method
---------------------
  BFS from root finds the first branch point (14082).  The two subtrees were
  characterised by their PVJ-node X centroids:
    Subtree → child 43979  (RBB): X_centroid = +0.051 m, 93 % on anatomic right
    Subtree → child 44015  (LBB): X_centroid = +0.011 m, 37 % on anatomic left
  The LBB fans across both the septum and the LV free wall — consistent with
  known anatomy and Strocchi biventricular model orientation.
"""

import re
import sys
import os

DRY_RUN = "--dry-run" in sys.argv

GRAPH_IN  = os.path.join(os.path.dirname(__file__), "constant", "purkinjeGraph")
GRAPH_OUT = os.path.join(os.path.dirname(__file__), "constant", "purkinjeGraph.lbbb")

# The single bridge edge that must be zeroed for complete LBBB.
LBB_BRIDGE_EDGE_INDEX = 44002          # 4(14082 44015 0.000251017 1)
LBB_BRIDGE_EXPECTED   = "4(14082 44015"  # sanity-check prefix

print(f"Reading  : {GRAPH_IN}")

with open(GRAPH_IN, "r") as f:
    content = f.read()

# ── locate every edge line in the conductionEdges block ─────────────────────
# Format: 4(nodeA nodeB length conductance)
edge_pattern = re.compile(r'4\(\d+ \d+ \S+ \S+\)')
edge_matches  = list(edge_pattern.finditer(content))

n_edges = len(edge_matches)
print(f"Edges found : {n_edges}")

if n_edges != 44050:
    print(f"WARNING: expected 44050 edges, found {n_edges}. Check the graph file.")

target = edge_matches[LBB_BRIDGE_EDGE_INDEX]
original_text = target.group(0)
print(f"Edge {LBB_BRIDGE_EDGE_INDEX}: {original_text}")

if not original_text.startswith(LBB_BRIDGE_EXPECTED):
    raise RuntimeError(
        f"Safety check failed: edge {LBB_BRIDGE_EDGE_INDEX} is '{original_text}', "
        f"expected it to start with '{LBB_BRIDGE_EXPECTED}'. "
        f"Re-run the identification script before proceeding."
    )

# ── build replacement: set 4th field (conductance) to 0 ─────────────────────
# 4(nodeA nodeB length conductance)  →  4(nodeA nodeB length 0)
m = re.match(r'(4\(\d+ \d+ \S+ )(\S+)(\))', original_text)
if not m:
    raise RuntimeError(f"Could not parse edge text: '{original_text}'")

replacement_text = m.group(1) + "0" + m.group(3)
print(f"Replacement : {replacement_text}")

if DRY_RUN:
    print("[dry-run] No file written.")
    sys.exit(0)

# Replace the single occurrence at the exact character position
new_content = (
    content[:target.start()]
    + replacement_text
    + content[target.end():]
)

with open(GRAPH_OUT, "w") as f:
    f.write(new_content)

print(f"Written  : {GRAPH_OUT}")
print()
print("To run the LBBB simulation:")
print("  cp constant/purkinjeGraph constant/purkinjeGraph.healthy")
print("  cp constant/purkinjeGraph.lbbb constant/purkinjeGraph")
print("  ./Allrun")
