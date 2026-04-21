import os
import re
import sys

def get_includes(file_path):
    includes = []
    # Standard C++ include regex: #include "file.H" or #include <file.H>
    include_pattern = re.compile(r'^\s*#include\s+["<]([^">]+)[">]')
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                match = include_pattern.match(line)
                if match:
                    includes.append(match.group(1))
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    
    return includes

def main(root_dir):
    if not os.path.exists(root_dir):
        print(f"Directory {root_dir} does not exist.")
        return

    print(f"Analyzing includes in: {root_dir}")
    print("-" * 50)

    for root, dirs, files in os.walk(root_dir):
        # Sort for consistent output
        dirs.sort()
        files.sort()
        
        # Skip some common build or hidden directories
        if 'Make' in dirs:
            dirs.remove('Make')
        if 'lnInclude' in dirs:
            dirs.remove('lnInclude')
            
        for file in files:
            if file.endswith(('.C', '.H', '.cpp', '.hpp')):
                file_path = os.path.join(root, file)
                rel_path = os.path.relpath(file_path, root_dir)
                includes = get_includes(file_path)
                
                if includes:
                    print(f"\nFile: {rel_path}")
                    for inc in includes:
                        print(f"  - {inc}")

if __name__ == "__main__":
    target_dir = "src/electroModels"
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    
    main(target_dir)
