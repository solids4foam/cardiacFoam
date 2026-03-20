from __future__ import annotations

import argparse
import json
from pathlib import Path

from .core.runtime.engine import DriverEngine
from .core.runtime.registry import list_tutorials, load_tutorial_spec
from .introspection import describe_tutorial
from .specs.common import default_setup_dir_name


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic OpenFOAM tutorial automation driver")
    parser.add_argument(
        "action",
        choices=["sim", "post", "all", "describe"],
        help="Pipeline stage to execute",
    )
    parser.add_argument(
        "--tutorial",
        required=True,
        help=(
            "Tutorial spec or case folder to run "
            f"({', '.join(list_tutorials())}, genericCase)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Plan and print simulation cases without running OpenFOAM.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue executing remaining cases after a failure.",
    )
    parser.add_argument(
        "--config",
        help=(
            "Path to JSON file with make_spec overrides. Supports either a top-level "
            "tutorial map (keys: singleCell, niederer2012, manufacturedFDA, "
            "restitutionCurves, genericCase/randomCase) or a direct parameter "
            "object for the selected tutorial."
        ),
    )
    parser.add_argument(
        "--tutorials-root",
        help=(
            "Optional path to the tutorials folder. Defaults to '<repo>/tutorials' when present."
        ),
    )
    return parser


def _load_spec_overrides(config_path: str, tutorial: str) -> dict:
    payload = json.loads(Path(config_path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("Config file must contain a JSON object")

    normalized_requested = tutorial.strip().casefold()
    for key, value in payload.items():
        if key.casefold() == normalized_requested:
            if not isinstance(value, dict):
                raise ValueError(f"Config section '{key}' must be a JSON object")
            return _normalize_spec_overrides(value)

    known_tutorial_keys = {
        *(name.casefold() for name in list_tutorials()),
        "genericcase",
        "randomcase",
    }
    if any(key.casefold() in known_tutorial_keys for key in payload):
        raise KeyError(
            f"No config section found for tutorial '{tutorial}'. "
            f"Available config sections: {', '.join(payload.keys())}"
        )

    return _normalize_spec_overrides(payload)


def _normalize_spec_overrides(overrides: dict) -> dict:
    normalized = dict(overrides)

    legacy_tutorial_name = normalized.pop("tutorial_name", None)
    if legacy_tutorial_name is not None:
        legacy_tutorial_name = str(legacy_tutorial_name)
        case_dir_name = normalized.get("case_dir_name")
        if case_dir_name is None:
            normalized["case_dir_name"] = legacy_tutorial_name
        elif str(case_dir_name) != legacy_tutorial_name:
            raise ValueError(
                "Config keys 'tutorial_name' and 'case_dir_name' must match. "
                "Use only 'case_dir_name'."
            )

    case_dir_name = normalized.get("case_dir_name")
    setup_dir_name = normalized.get("setup_dir_name")
    if case_dir_name is not None and setup_dir_name is not None:
        if str(setup_dir_name) == default_setup_dir_name(str(case_dir_name)):
            normalized.pop("setup_dir_name")

    return normalized


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.action == "post" and args.dry_run:
        parser.error("--dry-run is not valid with action=post")
    if args.action == "describe" and args.dry_run:
        parser.error("--dry-run is not valid with action=describe")
    if args.action == "describe" and args.continue_on_error:
        parser.error("--continue-on-error is not valid with action=describe")

    overrides = _load_spec_overrides(args.config, args.tutorial) if args.config else None
    if args.tutorials_root:
        if overrides is None:
            overrides = {}
        overrides["tutorials_root"] = args.tutorials_root

    if args.action == "describe":
        print(
            json.dumps(
                describe_tutorial(
                    args.tutorial,
                    overrides=overrides,
                    config_path=args.config,
                ),
                indent=2,
            )
        )
        return 0

    spec = load_tutorial_spec(args.tutorial, overrides=overrides)
    engine = DriverEngine(
        spec=spec,
        dry_run=args.dry_run,
        continue_on_error=args.continue_on_error,
        requested_action=args.action,
    )

    if args.action == "sim":
        engine.run_simulations()
    elif args.action == "post":
        engine.run_postprocess()
    else:
        engine.run_all()

    return 0
