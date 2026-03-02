from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from openfoam_driver.postprocessing.driver import PostprocessTask, run_postprocess_tasks


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

            self.assertEqual(manifest["schema_version"], "1.0")
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


if __name__ == "__main__":
    unittest.main()
