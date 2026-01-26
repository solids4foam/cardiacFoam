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
def parse_field_arrays(lines, i, current_section, drop_arrays=None):
    converted_blocks = []
    keep_arrays = []
    drop_arrays = drop_arrays or set()

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
        while i < len(lines) and lines[i].strip() == "":
            i += 1
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

        # Drop arrays explicitly requested
        if name in drop_arrays:
            print(f"Removed FIELD array '{name}' in section {current_section}.")
            continue

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
    current_section_count = None
    drop_scalar_arrays = {"vtkGhostType", "vtkOriginalPointIds", "vtkOriginalCellIds"}

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Track sections
        if stripped.startswith("POINT_DATA") or stripped.startswith("CELL_DATA"):
            tokens = stripped.split()
            current_section = tokens[0]
            try:
                current_section_count = int(tokens[1]) if len(tokens) > 1 else None
            except ValueError:
                current_section_count = None
            out.append(line)
            i += 1
            continue

        # Drop GLOBAL_IDS blocks (OpenFOAM reader does not support them)
        if stripped.startswith("GLOBAL_IDS"):
            if current_section_count is None:
                print("Warning: GLOBAL_IDS encountered without a valid POINT_DATA/CELL_DATA count; skipping line only.")
                i += 1
                continue
            i += 1
            _, i = read_values(lines, i, current_section_count)
            print(f"Removed GLOBAL_IDS block with {current_section_count} entries.")
            continue

        # Drop selected SCALARS arrays (e.g., vtkGhostType, vtkOriginalPointIds)
        if stripped.startswith("SCALARS"):
            tokens = stripped.split()
            name = tokens[1] if len(tokens) > 1 else ""
            if name in drop_scalar_arrays:
                i += 1  # skip SCALARS header
                if i < len(lines) and lines[i].strip().startswith("LOOKUP_TABLE"):
                    i += 1  # skip lookup table line
                if current_section_count is None:
                    print(f"Warning: SCALARS {name} encountered without a valid POINT_DATA/CELL_DATA count; skipping data read.")
                else:
                    _, i = read_values(lines, i, current_section_count)
                print(f"Removed SCALARS '{name}' with {current_section_count if current_section_count is not None else 'unknown'} entries.")
                continue

        # Handle FIELD
        if stripped.startswith("FIELD"):
            converted_blocks, keep_arrays, i = parse_field_arrays(
                lines,
                i,
                current_section,
                drop_arrays=drop_scalar_arrays,
            )
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
