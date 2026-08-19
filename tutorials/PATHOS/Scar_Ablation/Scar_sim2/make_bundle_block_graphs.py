#!/usr/bin/env python3
"""
make_bundle_block_graphs.py
===========================
Generates purkinjeGraph variants for bundle branch block simulations.

  constant/purkinjeGraph.lbbb  -- Left Bundle Branch Block  (edge 43954 zeroed)
  constant/purkinjeGraph.rbbb  -- Right Bundle Branch Block (edge 44002 zeroed)

Graph topology (His bundle bifurcation at node 14082):

  root (14083) --- His bundle (14082)
                        |
         edge 43954     |     edge 44002
              |                    |
        child 43979           child 44015
         LBB subtree            RBB subtree
         22158 nodes            21888 nodes

Usage
-----
  python3 make_bundle_block_graphs.py          # generate both files
  python3 make_bundle_block_graphs.py --lbbb   # generate + activate LBBB
  python3 make_bundle_block_graphs.py --rbbb   # generate + activate RBBB
  python3 make_bundle_block_graphs.py --dry-run
"""

import re, sys, os, shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GRAPH_IN   = os.path.join(SCRIPT_DIR, "constant", "purkinjeGraph")
OUT_LBBB   = os.path.join(SCRIPT_DIR, "constant", "purkinjeGraph.lbbb")
OUT_RBBB   = os.path.join(SCRIPT_DIR, "constant", "purkinjeGraph.rbbb")

DRY_RUN    = "--dry-run" in sys.argv
ACT_LBBB   = "--lbbb"    in sys.argv
ACT_RBBB   = "--rbbb"    in sys.argv

# ── bridge edges at the His bifurcation (node 14082) ────────────────────────
BRIDGES = {
    "lbbb": {
        "edge_index": 43954,
        "expected_prefix": "4(14082 43979",
        "out": OUT_LBBB,
        "description": "Left Bundle Branch Block  (LBB: child 43979)",
    },
    "rbbb": {
        "edge_index": 44002,
        "expected_prefix": "4(14082 44015",
        "out": OUT_RBBB,
        "description": "Right Bundle Branch Block (RBB: child 44015)",
    },
}

# Always build from the healthy base so modifications are never compounded.
# If a .healthy backup exists, use it; otherwise use the active graph.
GRAPH_BASE = GRAPH_IN + ".healthy" if os.path.exists(GRAPH_IN + ".healthy") else GRAPH_IN
print(f"Reading: {GRAPH_BASE}\n")
with open(GRAPH_BASE) as f:
    content = f.read()

edge_matches = list(re.compile(r'4\(\d+ \d+ \S+ \S+\)').finditer(content))
print(f"Total edges: {len(edge_matches)}")

for block_type, cfg in BRIDGES.items():
    idx     = cfg["edge_index"]
    prefix  = cfg["expected_prefix"]
    out     = cfg["out"]
    desc    = cfg["description"]

    target = edge_matches[idx]
    original = target.group(0)

    if not original.startswith(prefix):
        raise RuntimeError(
            f"Safety check failed for {block_type.upper()}: "
            f"edge {idx} is '{original}', expected prefix '{prefix}'."
        )

    m = re.match(r'(4\(\d+ \d+ \S+ )(\S+)(\))', original)
    replacement = m.group(1) + "0" + m.group(3)

    print(f"\n{desc}")
    print(f"  edge {idx}: {original}  -->  {replacement}")

    if DRY_RUN:
        print(f"  [dry-run] would write: {out}")
        continue

    new_content = content[:target.start()] + replacement + content[target.end():]
    with open(out, "w") as f:
        f.write(new_content)
    print(f"  Written: {out}")

if DRY_RUN:
    print("\n[dry-run] No files written.")
    sys.exit(0)

# ── optionally activate one variant ─────────────────────────────────────────
if ACT_LBBB:
    backup = GRAPH_IN + ".healthy"
    if not os.path.exists(backup):
        shutil.copy2(GRAPH_IN, backup)
        print(f"\nBacked up healthy graph to: {backup}")
    shutil.copy2(OUT_LBBB, GRAPH_IN)
    print(f"Activated LBBB: {OUT_LBBB} -> {GRAPH_IN}")
elif ACT_RBBB:
    backup = GRAPH_IN + ".healthy"
    if not os.path.exists(backup):
        shutil.copy2(GRAPH_IN, backup)
        print(f"\nBacked up healthy graph to: {backup}")
    shutil.copy2(OUT_RBBB, GRAPH_IN)
    print(f"Activated RBBB: {OUT_RBBB} -> {GRAPH_IN}")
else:
    print("\nNo variant activated. To activate:")
    print("  python3 make_bundle_block_graphs.py --lbbb")
    print("  python3 make_bundle_block_graphs.py --rbbb")
    print("To restore healthy:")
    print("  cp constant/purkinjeGraph.healthy constant/purkinjeGraph")
