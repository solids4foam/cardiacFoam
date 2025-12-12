from .vtk_utils import read_values



def write_wrapped(vals, per_line=9):
    out = []
    for j in range(0, len(vals), per_line):
        out.append(" ".join(vals[j:j+per_line]) + "\n")
    return out

def write_tensors(vals):
    out = []
    for t in range(0, len(vals), 9):
        a = vals[t:t+9]
        out.append(" ".join(a[0:3]) + "\n")
        out.append(" ".join(a[3:6]) + "\n")
        out.append(" ".join(a[6:9]) + "\n\n")
    return out

# --- Parse FIELD arrays inside a section ---
def parse_field_arrays(lines, i, current_section):
    converted_blocks = []
    keep_arrays = []

    field_header = lines[i].strip()
    tokens = field_header.split()
    if len(tokens) < 3:
        return converted_blocks, keep_arrays, i+1

    field_name = tokens[1]
    try:
        n_arrays = int(tokens[2])
    except ValueError:
        return converted_blocks, keep_arrays, i+1

    i += 1

    for arr_idx in range(n_arrays):
        if i >= len(lines):
            break
        hdr = lines[i].strip()
        parts = hdr.split()
        if len(parts) != 4:
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

        # Debug print
        print(f"\n--- FIELD ARRAY DEBUG ---")
        print(f"Current section: {current_section}, line:{i}")
        print(f"Array #{arr_idx+1}: {name}, nComp={ncomp}, nTuples={ntup}, type={dtype}")
        if current_section not in ("POINT_DATA", "CELL_DATA"):
            print(f"  WARNING: This array is outside POINT_DATA/CELL_DATA and will NOT be converted for OpenFOAM")

        # Convert if valid
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
        else:
            keep_arrays.append((f"{name} {ncomp} {ntup} {dtype}\n", vals))

    return converted_blocks, keep_arrays, i

# --- Main conversion function ---
def convert_vtk_file(infile, outfile):
    with open(infile, "r") as f:
        lines = f.readlines()

    out = []
    i = 0
    current_section = None

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Track sections
        if stripped.startswith("POINT_DATA") or stripped.startswith("CELL_DATA"):
            current_section = stripped.split()[0]
            out.append(line)
            i += 1
            continue

        # Handle FIELD
        if stripped.startswith("FIELD"):
            converted_blocks, keep_arrays, i = parse_field_arrays(lines, i, current_section)
            out.extend(converted_blocks)
            if keep_arrays:
                print(f"Arrays not converted in FIELD: {keep_arrays}")
                for hdr_line, vals in keep_arrays:
                    out.extend(write_wrapped(vals))
            continue

        # Default copy
        out.append(line)
        i += 1

    with open(outfile, "w") as f:
        f.writelines(out)

    print(f"Converted file written to {outfile}")

