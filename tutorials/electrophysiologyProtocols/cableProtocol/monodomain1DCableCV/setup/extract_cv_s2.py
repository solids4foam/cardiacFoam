#!/usr/bin/env python3

import sys
import re
from pathlib import Path

THRESHOLD = 0.0

def main():
    case_dir = Path.cwd()
    vm_files = list(case_dir.glob("postProcessing/cableProbes/*/Vm"))
    vm_files.sort(key=lambda x: float(x.parent.name))
    
    if not vm_files:
        print("FAILED: Vm files not found")
        sys.exit(1)
        
    positions = []
    probe_pattern = re.compile(r"#\s*Probe\s+\d+\s+\(([^)]+)\)")
    
    times = []
    vm_data = [] # List of lists: vm_data[probe_idx][time_idx]
    
    for vm_file in vm_files:
        with vm_file.open("r", encoding="ascii") as handle:
            for line in handle:
                line = line.strip()
                match = probe_pattern.match(line)
                if match:
                    coords = tuple(float(x) for x in match.group(1).split())
                    if len(positions) < 5:
                        positions.append(coords)
                    continue
                    
                if not line or line.startswith("#"):
                    continue
                    
                tokens = [float(x) for x in line.split()]
                if times and abs(tokens[0] - times[-1]) < 1e-8:
                    continue
                    
                times.append(tokens[0])
                
                if not vm_data:
                    vm_data = [[] for _ in range(len(tokens)-1)]
                    
                for i, val in enumerate(tokens[1:]):
                    vm_data[i].append(val)
                    
    num_probes = len(positions)
    
    crossings_per_probe = []
    for p_idx in range(num_probes):
        vms = vm_data[p_idx]
        crossings = []
        for t_idx in range(1, len(times)):
            # Only count crossings after t=4.30 to avoid plateau notches!
            if times[t_idx] < 4.30:
                continue
                
            if vms[t_idx-1] < THRESHOLD and vms[t_idx] >= THRESHOLD:
                v0 = vms[t_idx-1]
                v1 = vms[t_idx]
                t0 = times[t_idx-1]
                t1 = times[t_idx]
                t_cross = t0 + (THRESHOLD - v0) * (t1 - t0) / (v1 - v0)
                crossings.append(t_cross)
                
        clean_crossings = []
        if crossings:
            clean_crossings.append(crossings[0])
            for c in crossings[1:]:
                if c - clean_crossings[-1] > 0.050:
                    clean_crossings.append(c)
        crossings_per_probe.append(clean_crossings)
        
    if not crossings_per_probe[1] or not crossings_per_probe[3]:
        print("FAILED: S2 did not propagate")
        sys.exit(0)
        
    # S2 is just the very last crossing after 4.30
    t1_s2 = crossings_per_probe[1][-1]
    t3_s2 = crossings_per_probe[3][-1]
    
    dx = positions[3][0] - positions[1][0]
    dt = t3_s2 - t1_s2
    
    if dt <= 0:
        print("FAILED: dt <= 0")
        sys.exit(0)
        
    cv = dx / dt
    print(f"S2_CV = {cv:.4f} m/s")

if __name__ == "__main__":
    main()
