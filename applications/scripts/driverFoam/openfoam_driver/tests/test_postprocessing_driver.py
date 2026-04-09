from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks


def _repo_root_from_test() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "tutorials").exists() and (parent / "applications").exists():
            return parent
    raise RuntimeError("Could not locate repository root from test path")


class TestPostprocessingDriver(unittest.TestCase):
    def test_writes_schema_version_and_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            setup_root = root / "setup"
            output_dir = root / "out"
            setup_root.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)

            module_path = setup_root / "dummy_post.py"
            module_path.write_text(
                "\n".join(
                    [
                        "from pathlib import Path",
                        "def run_postprocessing(*, output_dir, setup_root=None, **kwargs):",
                        "    Path(output_dir).joinpath('plot_a.html').write_text('ok')",
                        "    return [",
                        "        {'path': 'plot_a.html', 'kind': 'plot', 'format': 'html', 'label': 'A'},",
                        "    ]",
                    ]
                )
            )

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="dummy",
                tasks=[PostprocessTask(module_relpath=Path("dummy_post.py"))],
            )

            manifest_path = output_dir / "plots.json"
            self.assertTrue(manifest_path.exists())
            manifest = json.loads(manifest_path.read_text())

            self.assertEqual(manifest["schema_version"], "1.1")
            self.assertEqual(manifest["artifact_count"], 1)
            artifact = manifest["artifacts"][0]
            self.assertEqual(artifact["path"], "plot_a.html")
            self.assertTrue(artifact["exists"])
            self.assertEqual(artifact["kind"], "plot")
            self.assertEqual(artifact["format"], "html")

    def test_strict_artifacts_raises_when_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            setup_root = root / "setup"
            output_dir = root / "out"
            setup_root.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)

            module_path = setup_root / "dummy_post.py"
            module_path.write_text(
                "\n".join(
                    [
                        "def run_postprocessing(*, output_dir, setup_root=None, **kwargs):",
                        "    return [{'path': 'missing_plot.html'}]",
                    ]
                )
            )

            with self.assertRaises(FileNotFoundError):
                run_postprocess_tasks(
                    setup_root=setup_root,
                    output_dir=output_dir,
                    tutorial_name="dummy",
                    strict_artifacts=True,
                    tasks=[PostprocessTask(module_relpath=Path("dummy_post.py"))],
                )

    def test_resolves_output_and_setup_placeholders(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            setup_root = root / "setup"
            output_dir = root / "out"
            setup_root.mkdir(parents=True, exist_ok=True)
            output_dir.mkdir(parents=True, exist_ok=True)

            module_path = setup_root / "dummy_post.py"
            module_path.write_text(
                "\n".join(
                    [
                        "import json",
                        "from pathlib import Path",
                        "def run_postprocessing(*, output_dir, setup_root=None, **kwargs):",
                        "    Path(output_dir).joinpath('kwargs.json').write_text(json.dumps(kwargs, sort_keys=True))",
                        "    Path(output_dir).joinpath('p.html').write_text('ok')",
                        "    return [{'path': 'p.html'}]",
                    ]
                )
            )

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="dummy",
                tasks=[
                    PostprocessTask(
                        module_relpath=Path("dummy_post.py"),
                        kwargs={
                            "output_token": "$OUTPUT_DIR/inside.txt",
                            "setup_token": "$SETUP_ROOT/here.txt",
                        },
                    )
                ],
            )

            kwargs_payload = json.loads((output_dir / "kwargs.json").read_text())
            self.assertEqual(kwargs_payload["output_token"], str(output_dir / "inside.txt"))
            self.assertEqual(kwargs_payload["setup_token"], str(setup_root / "here.txt"))

    def test_manufactured_postprocess_writes_csv_without_optional_plotting_deps(self) -> None:
        repo_root = _repo_root_from_test()
        setup_root = repo_root / "tutorials" / "manufacturedSolutions" / "monodomain" / "setupManufacturedFDA"

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            (output_dir / "1D_10_cells_implicit.dat").write_text(
                "\n".join(
                    [
                        "Vm 0 0 1.0e-02",
                        "u1 0 0 2.0e-02",
                        "u2 0 0 4.0e-02",
                    ]
                )
            )
            (output_dir / "1D_20_cells_implicit.dat").write_text(
                "\n".join(
                    [
                        "Vm 0 0 2.5e-03",
                        "u1 0 0 5.0e-03",
                        "u2 0 0 1.0e-02",
                    ]
                )
            )

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="manufacturedFDA",
                tasks=[PostprocessTask(module_relpath=Path("post_processing_manufactured.py"))],
            )

            csv_path = output_dir / "manufactured_convergence_rates.csv"
            manifest_path = output_dir / "plots.json"

            self.assertTrue(csv_path.exists())
            self.assertTrue(manifest_path.exists())
            self.assertIn("rate_Vm", csv_path.read_text())

            manifest = json.loads(manifest_path.read_text())
            artifact_paths = {artifact["path"] for artifact in manifest["artifacts"]}
            self.assertIn("manufactured_convergence_rates.csv", artifact_paths)

    def test_manufactured_postprocess_ignores_stale_outputs_not_in_manifest(self) -> None:
        repo_root = _repo_root_from_test()
        setup_root = repo_root / "tutorials" / "manufacturedSolutions" / "monodomain" / "setupManufacturedFDA"

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            (output_dir / "run_manifest.json").write_text(
                json.dumps(
                    {
                        "results": [
                            {
                                "status": "ok",
                                "params": {
                                    "dimension": "1D",
                                    "cells": 10,
                                    "solver": "implicit",
                                },
                            },
                            {
                                "status": "ok",
                                "params": {
                                    "dimension": "1D",
                                    "cells": 20,
                                    "solver": "implicit",
                                },
                            },
                        ]
                    }
                )
            )

            (output_dir / "1D_10_cells_implicit.dat").write_text(
                "\n".join(
                    [
                        "Vm 0 0 1.0e-02",
                        "u1 0 0 2.0e-02",
                        "u2 0 0 4.0e-02",
                    ]
                )
            )
            (output_dir / "1D_20_cells_implicit.dat").write_text(
                "\n".join(
                    [
                        "Vm 0 0 2.5e-03",
                        "u1 0 0 5.0e-03",
                        "u2 0 0 1.0e-02",
                    ]
                )
            )
            (output_dir / "1D_500_cells_implicit.dat").write_text(
                "\n".join(
                    [
                        "Vm 0 0 1.0e-08",
                        "u1 0 0 1.0e-08",
                        "u2 0 0 1.0e-08",
                    ]
                )
            )

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="manufacturedFDA",
                tasks=[PostprocessTask(module_relpath=Path("post_processing_manufactured.py"))],
            )

            csv_lines = (output_dir / "manufactured_convergence_rates.csv").read_text().strip().splitlines()
            self.assertEqual(len(csv_lines), 2)
            self.assertIn("10,20", csv_lines[1])


    def test_restitution_run_postprocessing_returns_list(self) -> None:
        """run_postprocessing must return list[dict], not None."""
        repo_root = _repo_root_from_test()
        module_path = (
            repo_root
            / "tutorials"
            / "restitutionCurves_s1s2Protocol"
            / "setupRestitutionCurves_s1s2Protocol"
            / "postProcessing_restCurves.py"
        )
        import importlib.util

        spec = importlib.util.spec_from_file_location("restcurves", module_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        fn = mod.run_postprocessing
        hints = fn.__annotations__
        return_hint = hints.get("return", None)
        # With `from __future__ import annotations` the hint is a string;
        # without it, it may be the actual type.  Accept both.
        is_list = return_hint is list or (
            isinstance(return_hint, str) and return_hint.startswith("list")
        )
        self.assertTrue(
            is_list,
            f"run_postprocessing must annotate return as list[dict], got: {return_hint!r}",
        )


    def test_niederer_table_summary_produces_csv_and_html(self) -> None:
        repo_root = _repo_root_from_test()
        setup_root = (
            repo_root / "tutorials" / "NiedererEtAl2012" / "setupNiedererEtAl2012"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            # Write a fake points CSV matching the expected pattern
            csv_content = "activationTime,Points:0,Points:1,Points:2\n0.042,0.0,0.0,0.007\n0.055,0.02,0.003,0.007\n"
            (output_dir / "implicit_TNNP_epicardialCells_points_DT001_DX01.csv").write_text(csv_content)

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="NiedererEtAl2012",
                tasks=[
                    PostprocessTask(
                        module_relpath=Path("postProcessing/table_summary.py")
                    )
                ],
            )

            self.assertTrue((output_dir / "NiedererEtAl2012_summary.csv").exists())
            self.assertTrue((output_dir / "NiedererEtAl2012_summary.html").exists())
            csv_text = (output_dir / "NiedererEtAl2012_summary.csv").read_text()
            self.assertIn("# tutorial: NiedererEtAl2012", csv_text)
            self.assertIn("case_id", csv_text)
            self.assertIn("implicit_TNNP_epicardialCells", csv_text)

    def test_singlecell_table_summary_produces_csv_and_html(self) -> None:
        repo_root = _repo_root_from_test()
        setup_root = repo_root / "tutorials" / "singleCell" / "setupSingleCell"

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            # Fake single-cell .txt output: time Vm (space-separated)
            # Resting ~-85 mV, peak ~40 mV, repolarises back to ~-85 mV
            import numpy as np
            t = np.linspace(0, 0.5, 500)
            vm = np.full_like(t, -85.0)
            vm[50:150] = np.linspace(-85, 40, 100)   # upstroke
            vm[150:350] = np.linspace(40, -85, 200)  # repolarisation
            txt_lines = ["time Vm"] + [f"{ti:.4f} {vi:.4f}" for ti, vi in zip(t, vm)]
            (output_dir / "TNNP_epicardialCells_run.txt").write_text("\n".join(txt_lines))

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="singleCell",
                tasks=[
                    PostprocessTask(
                        module_relpath=Path("postProcessing/table_summary.py")
                    )
                ],
            )

            self.assertTrue((output_dir / "singleCell_summary.csv").exists())
            self.assertTrue((output_dir / "singleCell_summary.html").exists())
            csv_text = (output_dir / "singleCell_summary.csv").read_text()
            self.assertIn("# tutorial: singleCell", csv_text)
            self.assertIn("APD_ms", csv_text)
            self.assertIn("peak_voltage_mV", csv_text)


    def test_restitution_table_summary_consolidates_model_csvs(self) -> None:
        repo_root = _repo_root_from_test()
        setup_root = (
            repo_root
            / "tutorials"
            / "restitutionCurves_s1s2Protocol"
            / "setupRestitutionCurves_s1s2Protocol"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            # Write fake per-model restitution CSVs
            (output_dir / "TNNP_restitution.csv").write_text(
                "tissue,DI_ms,APD90_ms\nepicardiaCells,300.0,280.0\nmCells,250.0,240.0\n"
            )
            (output_dir / "BuenoOrovio_restitution.csv").write_text(
                "tissue,DI_ms,APD90_ms\nepicardiaCells,310.0,290.0\n"
            )

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="restitutionCurves_s1s2Protocol",
                tasks=[
                    PostprocessTask(
                        module_relpath=Path("postProcessing/table_summary.py")
                    )
                ],
            )

            self.assertTrue((output_dir / "restitutionCurves_summary.csv").exists())
            self.assertTrue((output_dir / "restitutionCurves_summary.html").exists())
            csv_text = (output_dir / "restitutionCurves_summary.csv").read_text()
            self.assertIn("# tutorial: restitutionCurves_s1s2Protocol", csv_text)
            self.assertIn("ionic_model", csv_text)
            self.assertIn("TNNP", csv_text)
            self.assertIn("BuenoOrovio", csv_text)


if __name__ == "__main__":
    unittest.main()
