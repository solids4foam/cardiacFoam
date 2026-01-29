import shutil
from pathlib import Path
from string import Template


def instantiate_model_template(
    template_file: Path,
    model: str,
    year: int,
    outdir: Path,
    template_prefix: str,
):
    text = template_file.read_text(encoding="utf-8")

    # Support explicit placeholders and legacy tokens
    text = Template(text).safe_substitute(MODEL=model, YEAR=str(year))
    text = text.replace(template_prefix, model)
    text = text.replace("YYYY", str(year))

    new_name = template_file.name.replace(template_prefix, model)
    out_file = outdir / new_name
    out_file.write_text(text, encoding="utf-8")

    return out_file

def sort_folder(
    model: str,
    year: int,
    outdir: Path,
    model_kind: str = "ionic",
    verbose: bool = False,
):

    ionic_dir = outdir / f"{model}"
    ionic_dir.mkdir(parents=True, exist_ok=True)
    if verbose:
        print(f">>> Creating ionic folder: {ionic_dir}")

    script_dir = Path(__file__).resolve().parent
    template_dir = script_dir / "derivedClass_templates"
    template_prefix = "ionicGenericModel"
    if model_kind == "activeTension":
        template_prefix = "activeTensionGenericModel"

    templates = sorted(template_dir.glob(f"{template_prefix}.*"))
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
            template_prefix=template_prefix,
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
