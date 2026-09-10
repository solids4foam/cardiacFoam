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
#     sort_folder
#
# Description
#     Organizes generated code and instantiates model templates.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

import shutil
from pathlib import Path
from string import Template


def instantiate_model_template(
    template_file: Path,
    model: str,
    year: int,
    outdir: Path,
):
    text = template_file.read_text(encoding="utf-8")

    # Support explicit placeholders and legacy tokens
    text = Template(text).safe_substitute(MODEL=model, YEAR=str(year))
    text = text.replace("genericModel", model)
    text = text.replace("YYYY", str(year))

    new_name = template_file.name.replace("genericModel", model)
    out_file = outdir / new_name
    out_file.write_text(text, encoding="utf-8")

    return out_file

def sort_folder(
    model: str,
    year: int,
    outdir: Path,
    verbose: bool = False,
):

    ionic_dir = outdir / f"{model}"
    ionic_dir.mkdir(parents=True, exist_ok=True)
    if verbose:
        print(f">>> Creating ionic folder: {ionic_dir}")

    script_dir = Path(__file__).resolve().parent
    template_dir = script_dir / "derivedClass_templates"

    templates = sorted(template_dir.glob("genericModel.*"))
    if not templates:
        raise FileNotFoundError(f"No templates found in {template_dir}")

    created = []
    moved = []

    for tpl in templates:
        if verbose:
            print(f"    Instantiating template: {tpl.name}")

        created.append(instantiate_model_template(
            template_file=tpl,
            model=model,
            year=year,
            outdir=ionic_dir,
        ))

    # Copy generated headers
    generated_files = [
        outdir / f"{model}_{year}.H",
        outdir / f"{model}_{year}Names.H",
    ]

    for f in generated_files:
        if f.exists():
            if verbose:
                print(f"    Copying {f.name} → {ionic_dir}")
            shutil.move(f, ionic_dir / f.name)
            moved.append(ionic_dir / f.name)
        else:
            raise FileNotFoundError(f"Expected file not found: {f}")

    return {
        "ionic_dir": ionic_dir,
        "created": created,
        "moved": moved,
    }
