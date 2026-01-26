"""
convert_inspect_cli.py

Standalone CLI to inspect a VTK file and convert FIELD arrays to standard VTK data.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from utils.vtk_utils import inspect_fields, remove_blank_lines
from utils.vtk_convert_arrays_to_fields import convert_vtk_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect and convert VTK FIELD arrays.")
    parser.add_argument("--input", required=True, help="Input VTK file.")
    parser.add_argument("--output", default="output/converted.vtk", help="Output VTK file.")
    args = parser.parse_args()

    print(f"Inspecting input: {args.input}")
    inspect_fields(args.input)
    convert_vtk_file(args.input, args.output)
    remove_blank_lines(args.output, args.output)
    print(f"Converted file written to {args.output} (blank lines removed)")
    print(f"Inspecting output: {args.output}")
    inspect_fields(args.output)


if __name__ == "__main__":
    main()
