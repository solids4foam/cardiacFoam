import os
import re
from pathlib import Path

hierarchy = [
    'couplingModels',
    'genericWriter',
    'ionicModels',
    'activeTensionModels',
    'electroModels',
    'verificationModels',
    'electroMechanicalModels'
]

# Map each header basename to the libraries that provide it.
header_to_lib = {}
src_dir = Path(__file__).resolve().parent

for lib in hierarchy:
    lib_path = src_dir / lib
    if lib_path.exists():
        for root, _, files in os.walk(lib_path):
            if 'lnInclude' in root or 'Make' in root:
                continue
            for file in files:
                if file.endswith('.H'):
                    header_to_lib.setdefault(file, set()).add(lib)

violations = []

include_pattern = re.compile(r'^\s*#include\s+["<](.+?)[">]')

for lib_idx, lib in enumerate(hierarchy):
    lib_path = src_dir / lib
    if not lib_path.exists():
        continue

    illegal_libs = hierarchy[lib_idx+1:]

    for root, _, files in os.walk(lib_path):
        if 'lnInclude' in root or 'Make' in root:
            continue
        for file in files:
            if file.endswith('.H') or file.endswith('.C'):
                filepath = os.path.join(root, file)
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    for line_num, line in enumerate(f, 1):
                        match = include_pattern.match(line)
                        if match:
                            header_name = match.group(1).split('/')[-1]
                            if header_name in header_to_lib:
                                provider_libs = header_to_lib[header_name]
                                target_libs = set() if lib in provider_libs else provider_libs
                                illegal_targets = sorted(
                                    target_lib
                                    for target_lib in target_libs
                                    if target_lib in illegal_libs
                                )
                                if illegal_targets:
                                    violations.append({
                                        'file': os.path.relpath(filepath, src_dir),
                                        'line': line_num,
                                        'header': header_name,
                                        'from_lib': lib,
                                        'to_lib': ', '.join(illegal_targets)
                                    })

for v in violations:
    print(f"VIOLATION: {v['file']}:{v['line']} includes '{v['header']}' from higher-level library '{v['to_lib']}'")

if not violations:
    print("No hierarchical include violations found!")
