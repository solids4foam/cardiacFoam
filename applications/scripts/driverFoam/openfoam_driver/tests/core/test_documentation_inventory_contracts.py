"""Machine-checkable contracts for maintained model and utility inventories."""

import re
from pathlib import Path

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo


ROOT = Path(__file__).resolve().parents[6]


def _registered_types(library: str, table: str) -> set[str]:
    root = ROOT / "src" / library
    make_files = (root / "Make" / "files").read_text()
    sources = re.findall(r"^([^/#\s][^\s]*\.C)$", make_files, re.MULTILINE)
    registered: set[str] = set()

    pattern = re.compile(
        rf"addToRunTimeSelectionTable\s*\(\s*{table}\s*,\s*(\w+)\s*,",
        re.DOTALL,
    )
    for source in sources:
        text = (root / source).read_text()
        for class_name in pattern.findall(text):
            header = (root / source).with_suffix(".H").read_text()
            match = re.search(
                rf"class\s+{re.escape(class_name)}\b.*?"
                r"(?:OverrideTypeName|TypeName)\(\"([^\"]+)\"\)",
                header,
                re.DOTALL,
            )
            assert match, f"No runtime name found for {class_name} in {source}"
            registered.add(match.group(1))
    return registered


def _runtime_names_from_table(readme: Path) -> set[str]:
    text = readme.read_text()
    section = text.split("## Runtime model contract", 1)[1]
    section = section.split("\n## ", 1)[0]
    names = set()
    for line in section.splitlines():
        if line.startswith("|") and "`" in line:
            cells = [cell.strip() for cell in line.strip("|").split("|")]
            names.add(cells[1].strip("`"))
    return names


def test_ionic_runtime_inventory_matches_compiled_registrations():
    documented = _runtime_names_from_table(ROOT / "src/ionicModels/README.md")
    registered = _registered_types("ionicModels", "ionicModel")
    assert documented == registered
    assert all(
        not name.endswith("Batched") or name.endswith("compactBatched")
        for name in documented
    )


def test_active_tension_inventory_matches_compiled_registrations():
    documented = _runtime_names_from_table(
        ROOT / "src/activeTensionModels/README.md"
    )
    registered = _registered_types("activeTensionModels", "activeTensionModel")
    assert documented == registered
    assert "ManufacturedElectromechanics" in documented
    assert "NiedererHunterSmith" not in documented


def test_cuda_sources_cover_only_documented_batched_families():
    ionic_root = ROOT / "src/ionicModels"
    ionic_cuda = {
        Path(path).name.removesuffix("_cuda.cu")
        for path in re.findall(
            r"^([^#\s]+_cuda\.cu)$",
            (ionic_root / "Make/files-gpu").read_text(),
            re.MULTILINE,
        )
    }
    ionic_runtime = _runtime_names_from_table(ionic_root / "README.md")
    documented_batched_sources = {
        name.removesuffix("compactBatched") + "Batched"
        for name in ionic_runtime
        if name.endswith("compactBatched")
    }
    assert documented_batched_sources == ionic_cuda

    tension_root = ROOT / "src/activeTensionModels"
    tension_cuda = {
        Path(path).name.removesuffix("_cuda.cu")
        for path in re.findall(
            r"^([^#\s]+_cuda\.cu)$",
            (tension_root / "Make/files-gpu").read_text(),
            re.MULTILINE,
        )
    }
    tension_runtime = _runtime_names_from_table(tension_root / "README.md")
    assert {name for name in tension_runtime if name.endswith("Batched")} == tension_cuda


def test_utility_inventory_matches_make_targets_and_local_readmes():
    utilities = ROOT / "applications/utilities"
    expected = {}
    for make_file in utilities.glob("*/Make/files"):
        match = re.search(r"^EXE\s*=.*?/([^/\s]+)$", make_file.read_text(), re.MULTILINE)
        assert match, f"Missing EXE target in {make_file}"
        expected[make_file.parent.parent.name] = match.group(1)

    readme = (utilities / "README.md").read_text()
    rows = {}
    for line in readme.splitlines():
        if not line.startswith("|") or "[README](" not in line:
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        executable = cells[1].strip("`" )
        path_match = re.fullmatch(r"\[README\]\(([^/]+)/README\.md\)", cells[5])
        assert path_match, f"Malformed utility README link: {cells[5]}"
        rows[path_match.group(1)] = executable

    assert rows == expected
    for directory in rows:
        assert (utilities / directory / "README.md").is_file()
