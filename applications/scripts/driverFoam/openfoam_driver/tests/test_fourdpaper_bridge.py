from __future__ import annotations

import json
from pathlib import Path

from openfoam_driver.cli import main
from openfoam_driver.fourdpaper_bridge import export_4dpaper


def _write(path: Path, text: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _case(tmp_path: Path) -> Path:
    case = tmp_path / "tutorials" / "heartSim3D-1D" / "monodomainHeart"
    _write(case / "system" / "controlDict", "application cardiacFoam;\n")
    _write(case / "constant" / "electroProperties", "myocardiumSolver monodomainSolver;\n")
    _write(case / "monodomainHeart.foam", "")
    _write(case / "postProcessing" / "pseudoECG.dat", "# time V1\n0 0\n")
    _write(case / "postProcessing" / "pseudoECG_plots.png", "png")
    _write(
        case / "postProcessing" / "purkinjeNetworkVTK" / "purkinjeNetwork.vtk.series",
        '{"files": [{"name": "purkinjeNetwork_000001.vtk", "time": 0.0}]}',
    )
    return case


def _paper(tmp_path: Path) -> Path:
    paper = tmp_path / "paper"
    _write(paper / "main.qmd", "---\ntitle: Paper\nformat: html\n---\n")
    return paper


def test_export_4dpaper_stages_case_and_writes_qmd_manifest(tmp_path: Path):
    case = _case(tmp_path)
    paper = _paper(tmp_path)

    result = export_4dpaper(case_root=case, paper_root=paper, copy=True)

    assert result.status == "staged"
    assert (paper / "data" / "monodomainHeart" / "monodomainHeart.foam").is_file()
    qmd = Path(result.qmd_fragment).read_text()
    assert "4d-image" in qmd
    assert "4d-timeseries" in qmd
    assert "4d-panel" in qmd
    assert "4d-multi-image" in qmd
    assert "pseudoECG.dat" in qmd
    manifest = json.loads(Path(result.manifest_path).read_text())
    assert manifest["status"] == "staged"
    assert manifest["fields"] == ["Vm", "activationTime"]
    assert "sections/monodomainHeart_results.qmd" in (paper / "main.qmd").read_text()


def test_export_4dpaper_collects_rendered_assets(tmp_path: Path, monkeypatch):
    case = _case(tmp_path)
    paper = _paper(tmp_path)
    _write(paper / "state" / "figures" / "monodomainHeart-vm.png", "png")

    def fake_run(*_, **__):
        return ["docker", "compose"], 0, "ok", ""

    monkeypatch.setattr("openfoam_driver.fourdpaper_bridge._run_render", fake_run)

    result = export_4dpaper(case_root=case, paper_root=paper, render=True)

    assert result.status == "rendered"
    assert result.rendered_assets
    copied = case / result.rendered_assets[0].path
    assert copied.is_file()


def test_export_4dpaper_cli_prints_json(tmp_path: Path, capsys):
    case = _case(tmp_path)
    paper = _paper(tmp_path)

    code = main([
        "export-4dpaper",
        "--case", str(case),
        "--paper", str(paper),
        "--copy",
        "--no-update-main",
    ])

    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "staged"
    assert payload["data_path"] == "data/monodomainHeart"

