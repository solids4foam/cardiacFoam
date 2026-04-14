#!/usr/bin/env python3
"""
generate_catalog.py — CI tool to verify the static ionic_model_catalog.py
matches the runtime output of listCellModelsVariables.

Usage:
    python3 generate_catalog.py --report postProcessing/listCellModelsVariables.txt --model TNNP

Run after recompiling cardiacFoam to detect drift between C++ model variables
and the static Python catalog.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys
from pathlib import Path


def parse_report_text(text: str) -> dict[str, list[str]]:
    """Parse listCellModelsVariables.txt report from string.

    Returns dict with keys: 'ionic_model', 'states', 'algebraic', 'constants'
    Each value is a list of variable name strings.
    """
    result = {
        "ionic_model": None,
        "states": [],
        "algebraic": [],
        "constants": [],
    }

    # Extract ionic model name
    model_match = re.search(r"Selected ionicModel:\s*(\w+)", text)
    if model_match:
        result["ionic_model"] = model_match.group(1)

    # Parse constants section
    constants_section = re.search(
        r"Ionic constants\s*\(\d+\).*?\n((?:\s+constants\s*\[\d+\].*\n)*)",
        text,
        re.DOTALL
    )
    if constants_section:
        const_lines = constants_section.group(1)
        for line in const_lines.split("\n"):
            match = re.search(r"constants\s*\[\d+\]\s+(\w+)\s*-->", line)
            if match:
                result["constants"].append(match.group(1))

    # Parse states section
    states_section = re.search(
        r"Ionic states\s*\(\d+\).*?\n((?:\s+states\s*\[\d+\].*\n)*)",
        text,
        re.DOTALL
    )
    if states_section:
        state_lines = states_section.group(1)
        for line in state_lines.split("\n"):
            match = re.search(r"states\s*\[\d+\]\s+(\w+)\s*-->", line)
            if match:
                result["states"].append(match.group(1))

    # Parse algebraic section
    algebraic_section = re.search(
        r"Ionic algebraic\s*\(\d+\).*?\n((?:\s+algebraic\s*\[\d+\].*\n)*)",
        text,
        re.DOTALL
    )
    if algebraic_section:
        alg_lines = algebraic_section.group(1)
        for line in alg_lines.split("\n"):
            match = re.search(r"algebraic\s*\[\d+\]\s+(\w+)", line)
            if match:
                result["algebraic"].append(match.group(1))

    return result


def parse_report(report_path: Path) -> dict[str, list[str]]:
    """Parse a listCellModelsVariables.txt report file.

    Returns dict with keys: 'ionic_model', 'states', 'algebraic', 'constants'
    Each value is a list of variable name strings.
    """
    with open(report_path, "r") as f:
        text = f.read()
    return parse_report_text(text)


def compare_with_catalog(model_name: str, parsed: dict[str, list[str]]) -> list[str]:
    """Compare parsed report against static catalog.

    Returns list of difference strings (empty = no diff).
    """
    # Import the catalog module directly via sys.path manipulation
    catalog_dir = pathlib.Path(__file__).parent.parent
    if str(catalog_dir) not in sys.path:
        sys.path.insert(0, str(catalog_dir))

    from ionic_model_catalog import IONIC_MODEL_CATALOG
    catalog = IONIC_MODEL_CATALOG

    if model_name not in catalog:
        return [f"ERROR: model {model_name!r} not found in catalog"]

    entry = catalog[model_name]

    diffs = []

    # Compare states
    runtime_states = set(parsed["states"])
    catalog_states = set(entry.states)

    missing_in_catalog = runtime_states - catalog_states
    missing_in_runtime = catalog_states - runtime_states

    for var in sorted(missing_in_catalog):
        diffs.append(f"+ states: {var} (in runtime, not in catalog)")

    for var in sorted(missing_in_runtime):
        diffs.append(f"- states: {var} (in catalog, not in runtime)")

    # Compare algebraic
    runtime_alg = set(parsed["algebraic"])
    catalog_alg = set(entry.algebraic)

    missing_in_catalog = runtime_alg - catalog_alg
    missing_in_runtime = catalog_alg - runtime_alg

    for var in sorted(missing_in_catalog):
        diffs.append(f"+ algebraic: {var} (in runtime, not in catalog)")

    for var in sorted(missing_in_runtime):
        diffs.append(f"- algebraic: {var} (in catalog, not in runtime)")

    # Compare constants
    runtime_consts = set(parsed["constants"])
    catalog_consts = set(entry.constants)

    missing_in_catalog = runtime_consts - catalog_consts
    missing_in_runtime = catalog_consts - runtime_consts

    for var in sorted(missing_in_catalog):
        diffs.append(f"+ constants: {var} (in runtime, not in catalog)")

    for var in sorted(missing_in_runtime):
        diffs.append(f"- constants: {var} (in catalog, not in runtime)")

    return diffs


def print_update_snippet(model_name: str, parsed: dict[str, list[str]]) -> None:
    """Print a Python snippet showing what the catalog entry should look like."""
    states_str = ", ".join(f'"{s}"' for s in parsed["states"])
    algebraic_str = ", ".join(f'"{a}"' for a in parsed["algebraic"])
    constants_str = ", ".join(f'"{c}"' for c in parsed["constants"])

    print(f'    "{model_name}": IonicModelEntry(')
    print(f'        states=({states_str}),')
    print(f'        algebraic=({algebraic_str}),')
    print(f'        constants=({constants_str}),')
    print(f'        recommended_exports=(...),')
    print(f'        compatible_tissues=(...),')
    print(f'        compatible_solvers=(...),')
    print(f'        species=(...),')
    print(f'        cardiac_region=(...),')
    print(f'        model_type="ionic",')
    print(f'        description="...",')
    print(f'    ),')


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare listCellModelsVariables output with static catalog."
    )
    parser.add_argument(
        "--report",
        required=True,
        type=Path,
        help="Path to listCellModelsVariables.txt",
    )
    parser.add_argument(
        "--model",
        required=True,
        help="Ionic model name (e.g. TNNP)",
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="Print catalog update snippet",
    )
    args = parser.parse_args()

    if not args.report.exists():
        print(f"ERROR: report file not found: {args.report}", file=sys.stderr)
        return 1

    parsed = parse_report(args.report)

    # Validate the report matches the requested model
    if parsed.get("ionic_model") and parsed["ionic_model"] != args.model:
        print(
            f"WARNING: report says model={parsed['ionic_model']!r} but --model={args.model!r}",
            file=sys.stderr,
        )

    diffs = compare_with_catalog(args.model, parsed)

    if diffs:
        print(f"DIFF for {args.model}:")
        for d in diffs:
            print(f"  {d}")
        if args.update:
            print()
            print_update_snippet(args.model, parsed)
        return 1

    print(
        f"OK: {args.model} catalog matches runtime report "
        f"({len(parsed['states'])} states, {len(parsed['algebraic'])} algebraic, "
        f"{len(parsed['constants'])} constants)"
    )
    if args.update:
        print()
        print_update_snippet(args.model, parsed)
    return 0


# Quick self-test
def _self_test() -> None:
    """Verify parse_report_text correctly extracts variable names."""
    sample = """
========== listCellModelsVariables ==========

Selected ionicModel: TestModel

Ionic constants (2) --> initial value
  constants [0] AC_R --> 8314.0
  constants [1] AC_T --> 310.0

Ionic states (3) --> initial value
  states [0] Vm --> -85.23
  states [1] m --> 0.001
  states [2] h --> 0.9

Ionic algebraic (1)
  algebraic [0] Iion
"""
    result = parse_report_text(sample)
    assert result["ionic_model"] == "TestModel", f"Expected TestModel, got {result['ionic_model']}"
    assert result["states"] == ["Vm", "m", "h"], f"Expected ['Vm', 'm', 'h'], got {result['states']}"
    assert result["algebraic"] == ["Iion"], f"Expected ['Iion'], got {result['algebraic']}"
    assert result["constants"] == ["AC_R", "AC_T"], f"Expected ['AC_R', 'AC_T'], got {result['constants']}"
    print("Self-test PASSED")


if __name__ == "__main__":
    sys.exit(main())
