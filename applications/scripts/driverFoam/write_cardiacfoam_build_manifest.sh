#!/usr/bin/env bash
set -euo pipefail

# Write the build identity consumed by the cardiacFoam driverFOAM plugin.
# The script is called by the repository's top-level Allwmake after all
# libraries and applications have compiled successfully.

if [[ -z "${WM_PROJECT_DIR:-}" || -z "${WM_OPTIONS:-}" ]]; then
    echo "OpenFOAM must be sourced before writing the cardiacFoam build manifest." >&2
    exit 1
fi

manifest_dir="${FOAM_USER_LIBBIN:-${FOAM_LIBBIN:-}}"
if [[ -z "$manifest_dir" ]]; then
    echo "FOAM_USER_LIBBIN or FOAM_LIBBIN is required." >&2
    exit 1
fi
mkdir -p "$manifest_dir"
manifest_path="$manifest_dir/cardiacFoam.build.json"

backend="${DRIVERFOAM_CARDIACFOAM_BACKEND:-}"
if [[ "$backend" != "full" && "$backend" != "lightweight" ]]; then
    echo "DRIVERFOAM_CARDIACFOAM_BACKEND must be 'full' or 'lightweight'." >&2
    exit 1
fi

solids_root="${DRIVERFOAM_CARDIACFOAM_SOLIDS4FOAM_ROOT:-${SOLIDS4FOAM_INST_DIR:-}}"
if [[ "$backend" == "full" && -z "$solids_root" ]]; then
    echo "Full backend requires DRIVERFOAM_CARDIACFOAM_SOLIDS4FOAM_ROOT." >&2
    exit 1
fi

user_appbin="${FOAM_USER_APPBIN:-}"
solver="$user_appbin/cardiacFoam"
if [[ ! -f "$solver" ]]; then
    echo "Compiled cardiacFoam not found at $solver." >&2
    exit 1
fi

library_dirs=()
for candidate in "${FOAM_USER_LIBBIN:-}" "${FOAM_MODULE_LIBBIN:-}" "${FOAM_LIBBIN:-}"; do
    [[ -n "$candidate" ]] && library_dirs+=("$candidate")
done

find_library() {
    local name="$1"
    local directory
    for directory in "${library_dirs[@]}"; do
        for extension in dylib so; do
            if [[ -f "$directory/lib${name}.${extension}" ]]; then
                printf '%s\n' "$directory/lib${name}.${extension}"
                return 0
            fi
        done
    done
    return 1
}

required_libraries=(electroModels ionicModels genericWriter activeTensionModels)
if [[ "$backend" == "full" ]]; then
    required_libraries+=(solids4FoamModels electroMechanicalModels)
else
    required_libraries+=(physicsModel)
fi

artifact_paths=("$solver")
artifact_names=("cardiacFoam")
for library in "${required_libraries[@]}"; do
    path="$(find_library "$library" || true)"
    if [[ -z "$path" ]]; then
        echo "Required library lib${library} was not found." >&2
        exit 1
    fi
    artifact_paths+=("$path")
    artifact_names+=("lib${library}")
done

linked_libraries=()
if command -v otool >/dev/null 2>&1; then
    while IFS= read -r library; do
        [[ -n "$library" ]] && linked_libraries+=("$library")
    done < <(otool -L "$solver" | awk 'NR > 1 {print $1}')
fi

source_revision=""
if [[ -n "$solids_root" && -e "$solids_root/.git" ]]; then
    source_revision="$(git -C "$solids_root" rev-parse HEAD 2>/dev/null || true)"
fi

python3 - "$manifest_path" "$backend" "$WM_PROJECT_DIR" "${WM_PROJECT_VERSION:-}" "$WM_OPTIONS" "$solids_root" "$source_revision" "${artifact_names[*]}" "${artifact_paths[*]}" "${linked_libraries[*]}" <<'PY'
import hashlib
import json
import pathlib
import sys

manifest_path, backend, foam_root, foam_version, wm_options, solids_root, source_revision, names, paths, linked = sys.argv[1:]
name_list = names.split()
path_list = paths.split()
linked_list = linked.split()
artifacts = []
for name, path in zip(name_list, path_list):
    data = pathlib.Path(path).read_bytes()
    artifacts.append({"name": name, "path": path, "sha256": hashlib.sha256(data).hexdigest()})
payload = {
    "schema_version": 1,
    "plugin": "org.cardiacfoam",
    "backend": backend,
    "openfoam": {"root": str(pathlib.Path(foam_root).resolve()), "version": foam_version, "options": wm_options},
    "solids4foam": {"root": str(pathlib.Path(solids_root).resolve()) if solids_root else None, "revision": source_revision or None},
    "linked_libraries": linked_list,
    "artifacts": artifacts,
}
pathlib.Path(manifest_path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY

echo "Wrote cardiacFoam build manifest: $manifest_path"
