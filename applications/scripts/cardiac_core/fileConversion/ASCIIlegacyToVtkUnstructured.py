"""
convert_inspect_cli.py

Standalone CLI to inspect a VTK file and convert FIELD arrays to standard VTK data.
"""
import argparse

from cardiac_core.io.postprocess import inspect_vtk, postprocess_vtk_output
from cardiac_core.utils.vtk_convert_arrays_to_fields import convert_vtk_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect and convert VTK FIELD arrays.")
    parser.add_argument("--input", required=True, help="Input VTK file.")
    parser.add_argument("--output", default="outputs/converted.vtk", help="Output VTK file.")
    args = parser.parse_args()

    print(f"Inspecting input: {args.input}")
    inspect_vtk(args.input)
    convert_vtk_file(args.input, args.output)
    postprocess_vtk_output(args.output, convert_fields=False, remove_blanks=True)
    print(f"Converted file written to {args.output} (blank lines removed)")
    print(f"Inspecting output: {args.output}")
    inspect_vtk(args.output)


if __name__ == "__main__":
    main()
