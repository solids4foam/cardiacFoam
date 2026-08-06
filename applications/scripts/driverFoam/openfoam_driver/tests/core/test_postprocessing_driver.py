#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     test_postprocessing_driver
#
# Description
#     Tests postprocessing driver logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import json
import importlib.util
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks
from openfoam_driver.tests.conftest import monorepo_root, skip_without_monorepo


@skip_without_monorepo

class TestPostprocessingDriver(unittest.TestCase):
    def test_monodomain_1d_cable_cv_postprocess_writes_csv_and_plots(self) -> None:
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "electrophysiologyProtocols/cableProtocol"
            / "monodomain1DCableCV"
            / "setup"
        )

        module_path = setup_root / "table_summary.py"
        spec = importlib.util.spec_from_file_location("monodomain1d_table_summary", module_path)
        if spec is None or spec.loader is None:
            raise RuntimeError("Could not load monodomain1DCableCV table_summary module")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        has_matplotlib = bool(module._has_matplotlib())

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            cases = [
                ("implicit_BuenoOrovio_epicardialCells_DT0.005_DX0.1_COND01", 0.6951652, 0.01438507),
                ("implicit_BuenoOrovio_epicardialCells_DT0.01_DX0.1_COND01", 0.68399827, 0.01461992),
                ("implicit_BuenoOrovio_epicardialCells_DT0.005_DX0.2_COND01", 0.65291796, 0.01531586),
                ("implicit_BuenoOrovio_epicardialCells_DT0.01_DX0.2_COND01", 0.6433723, 0.01554310),
            ]
            for case_id, cv_value, dt_s in cases:
                payload = {
                    "case_id": case_id,
                    "central_cv": {
                        "start_probe": 1,
                        "end_probe": 3,
                        "dx_m": 0.01,
                        "dt_s": dt_s,
                        "cv_m_per_s": cv_value,
                    },
                }
                (output_dir / f"{case_id}_cv_summary.json").write_text(json.dumps(payload))

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="monodomainAndEikonal1DCableCVConvergence",
                tasks=[PostprocessTask(module_relpath=Path("table_summary.py"))],
            )

            csv_path = output_dir / "monodomainAndEikonal1DCableCVConvergence_summary.csv"
            html_path = output_dir / "monodomainAndEikonal1DCableCVConvergence_summary.html"
            manifest = json.loads((output_dir / "plots.json").read_text())

            self.assertTrue(csv_path.exists())
            self.assertFalse(html_path.exists())

            artifact_paths = {artifact["path"] for artifact in manifest["artifacts"]}
            self.assertIn("monodomainAndEikonal1DCableCVConvergence_summary.csv", artifact_paths)
            self.assertNotIn("monodomainAndEikonal1DCableCVConvergence_summary.html", artifact_paths)

            plot_artifacts = [artifact for artifact in manifest["artifacts"] if artifact["kind"] == "plot"]
            if has_matplotlib:
                self.assertEqual(len(plot_artifacts), 2)
                for artifact in plot_artifacts:
                    self.assertEqual(artifact["format"], "png")
                    self.assertTrue((output_dir / artifact["path"]).exists())
            else:
                self.assertEqual(len(plot_artifacts), 0)

    def test_monodomain_1d_cable_cv_postprocess_supports_ionic_model_subfolders(self) -> None:
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "electrophysiologyProtocols/cableProtocol"
            / "monodomain1DCableCV"
            / "setup"
        )

        module_path = setup_root / "table_summary.py"
        spec = importlib.util.spec_from_file_location("monodomain1d_table_summary_nested", module_path)
        if spec is None or spec.loader is None:
            raise RuntimeError("Could not load monodomain1DCableCV table_summary module")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        has_matplotlib = bool(module._has_matplotlib())

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            model_dir = output_dir / "BuenoOrovio"
            model_dir.mkdir(parents=True, exist_ok=True)

            cases = [
                ("implicit_BuenoOrovio_epicardialCells_DT0.005_DX0.1_COND01", 0.6951652, 0.01438507),
                ("implicit_BuenoOrovio_epicardialCells_DT0.01_DX0.1_COND01", 0.68399827, 0.01461992),
                ("implicit_BuenoOrovio_epicardialCells_DT0.005_DX0.2_COND01", 0.65291796, 0.01531586),
                ("implicit_BuenoOrovio_epicardialCells_DT0.01_DX0.2_COND01", 0.6433723, 0.01554310),
            ]
            for case_id, cv_value, dt_s in cases:
                payload = {
                    "case_id": case_id,
                    "central_cv": {
                        "start_probe": 1,
                        "end_probe": 3,
                        "dx_m": 0.01,
                        "dt_s": dt_s,
                        "cv_m_per_s": cv_value,
                    },
                }
                (model_dir / f"{case_id}_cv_summary.json").write_text(json.dumps(payload))

            run_postprocess_tasks(
                setup_root=setup_root,
                output_dir=output_dir,
                tutorial_name="monodomainAndEikonal1DCableCVConvergence",
                tasks=[PostprocessTask(module_relpath=Path("table_summary.py"))],
            )

            csv_path = model_dir / "monodomainAndEikonal1DCableCVConvergence_summary.csv"
            html_path = model_dir / "monodomainAndEikonal1DCableCVConvergence_summary.html"
            manifest = json.loads((output_dir / "plots.json").read_text())

            self.assertTrue(csv_path.exists())
            self.assertFalse(html_path.exists())

            artifact_paths = {artifact["path"] for artifact in manifest["artifacts"]}
            self.assertIn("BuenoOrovio/monodomainAndEikonal1DCableCVConvergence_summary.csv", artifact_paths)
            self.assertNotIn("monodomainAndEikonal1DCableCVConvergence_summary.csv", artifact_paths)

            plot_artifacts = [artifact for artifact in manifest["artifacts"] if artifact["kind"] == "plot"]
            if has_matplotlib:
                self.assertEqual(len(plot_artifacts), 2)
                for artifact in plot_artifacts:
                    self.assertTrue(str(artifact["path"]).startswith("BuenoOrovio/"))
                    self.assertEqual(artifact["format"], "png")
                    self.assertTrue((output_dir / artifact["path"]).exists())
            else:
                self.assertEqual(len(plot_artifacts), 0)

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
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "manufacturedSolutions"
            / "monodomainPseudoECG"
            / "setup"
        )

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
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "manufacturedSolutions"
            / "monodomainPseudoECG"
            / "setup"
        )

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
        repo_root = monorepo_root  # type: ignore[assignment]
        module_path = (
            repo_root
            / "tutorials"
            / "electrophysiologyProtocols"
            / "restitutionCurves_s1s2Protocol"
            / "setup"
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


    def test_singlecell_table_summary_produces_csv_and_html(self) -> None:
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "electrophysiologyProtocols"
            / "singleCell"
            / "setup"
        )

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
                        module_relpath=Path("table_summary.py")
                    )
                ],
            )

            self.assertTrue((output_dir / "singleCell_summary.csv").exists())
            self.assertTrue((output_dir / "singleCell_summary.html").exists())
            csv_text = (output_dir / "singleCell_summary.csv").read_text()
            self.assertIn("# entry: singleCell", csv_text)
            self.assertIn("APD_ms", csv_text)
            self.assertIn("peak_voltage_mV", csv_text)


    def test_restitution_table_summary_consolidates_model_csvs(self) -> None:
        repo_root = monorepo_root  # type: ignore[assignment]
        setup_root = (
            repo_root
            / "tutorials"
            / "electrophysiologyProtocols"
            / "restitutionCurves_s1s2Protocol"
            / "setup"
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
                        module_relpath=Path("table_summary.py")
                    )
                ],
            )

            self.assertTrue((output_dir / "restitutionCurves_summary.csv").exists())
            self.assertTrue((output_dir / "restitutionCurves_summary.html").exists())
            csv_text = (output_dir / "restitutionCurves_summary.csv").read_text()
            self.assertIn("# entry: restitutionCurves_s1s2Protocol", csv_text)
            self.assertIn("ionic_model", csv_text)
            self.assertIn("TNNP", csv_text)
            self.assertIn("BuenoOrovio", csv_text)


if __name__ == "__main__":
    unittest.main()
