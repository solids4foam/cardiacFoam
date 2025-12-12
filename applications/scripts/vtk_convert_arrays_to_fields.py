#!/usr/bin/env python3
# Convert legacy VTK FIELD arrays (1/3/9 comps) to SCALARS/VECTORS/TENSORS.

infile = "ASCIIlegacy.vtk"
outfile = "ASCIIlegacy_fixed.vtk"

def read_values(lines, i, nvals):
    vals = []
    while i < len(lines) and len(vals) < nvals:
        vals.extend(lines[i].split())
        i += 1
    if len(vals) != nvals:
        raise ValueError(f"Expected {nvals} values, got {len(vals)}")
    return vals, i

def write_wrapped(vals, per_line=9):
    out = []
    for j in range(0, len(vals), per_line):
        out.append(" ".join(vals[j:j+per_line]) + "\n")
    return out

def write_tensors(vals):
    # vals length is 9 * ntuples; write each tensor as 3 rows of 3
    out = []
    for t in range(0, len(vals), 9):
        a = vals[t:t+9]
        out.append(" ".join(a[0:3]) + "\n")
        out.append(" ".join(a[3:6]) + "\n")
        out.append(" ".join(a[6:9]) + "\n\n")
    return out

with open(infile, "r") as f:
    lines = f.readlines()

out = []
i = 0
current_section = None  # None | "POINT_DATA" | "CELL_DATA"

while i < len(lines):
    line = lines[i]
    stripped = line.strip()

    # Track whether we're under POINT_DATA or CELL_DATA (needed to place attributes)
    if stripped.startswith("POINT_DATA") or stripped.startswith("CELL_DATA"):
        current_section = stripped.split()[0]  # "POINT_DATA" or "CELL_DATA"
        out.append(line)
        i += 1
        continue

    # Handle FIELD blocks
    if stripped.startswith("FIELD"):
        # Example: FIELD FieldData 3
        tokens = stripped.split()
        # Some files use: FIELD FieldData <nArrays>
        # Keep original header tokens around, but we will rebuild if needed.
        if len(tokens) < 3:
            out.append(line); i += 1; continue
        field_name = tokens[1]
        try:
            n_arrays = int(tokens[2])
        except ValueError:
            out.append(line); i += 1; continue

        i += 1

        # Parse the following n_arrays arrays
        keep_arrays = []     # list of (header_line, values_as_list)
        converted_blocks = []  # list of strings for converted arrays to emit

        for _ in range(n_arrays):
            hdr = lines[i].strip()
            parts = hdr.split()
            # Expect: <name> <nComp> <nTuples> <type>
            if len(parts) != 4:
                # Malformed; keep as-is and move on one line
                keep_arrays.append((lines[i], []))
                i += 1
                continue

            name, ncomp_s, ntup_s, dtype = parts
            try:
                ncomp = int(ncomp_s)
                ntup = int(ntup_s)
            except ValueError:
                keep_arrays.append((lines[i], []))
                i += 1
                continue

            i += 1
            nvals = ncomp * ntup
            vals, i = read_values(lines, i, nvals)

            # We only convert if we're under POINT_DATA or CELL_DATA (valid locations)
            if current_section in ("POINT_DATA", "CELL_DATA") and ncomp in (1, 3, 9):
                if ncomp == 1:
                    converted_blocks.append(f"SCALARS {name} {dtype}\n")
                    converted_blocks.append("LOOKUP_TABLE default\n")
                    converted_blocks.extend(write_wrapped(vals))
                elif ncomp == 3:
                    converted_blocks.append(f"VECTORS {name} {dtype}\n")
                    converted_blocks.extend(write_wrapped(vals))
                elif ncomp == 9:
                    converted_blocks.append(f"TENSORS {name} {dtype}\n")
                    converted_blocks.extend(write_tensors(vals))
                # Do NOT keep this one in FIELD; it's converted.
            else:
                # Keep untouched inside FIELD
                keep_arrays.append((f"{name} {ncomp} {ntup} {dtype}\n", vals))

        # Emit converted arrays first (in place of the FIELD block)
        out.extend(converted_blocks)

        # Rebuild the FIELD block if any arrays remain
        if keep_arrays:
            out.append(f"FIELD {field_name} {len(keep_arrays)}\n")
            for hdr_line, vals in keep_arrays:
                out.append(hdr_line if hdr_line.endswith("\n") else hdr_line + "\n")
                if vals:
                    out.extend(write_wrapped(vals))
        # If none remain, we drop the FIELD block entirely.
        continue

    # Default: copy line through
    out.append(line)
    i += 1

with open(outfile, "w") as f:
    f.writelines(out)

print(f"Converted file written to {outfile}")
