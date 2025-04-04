import json
import shutil
from pathlib import Path

import numpy as np

# import cardiac_benchmark
# import matplotlib.pyplot as plt

step2dir = Path("/global/D1/homes/henriknf/cardiac_benchmark/step2")


TABLE = r"""
\begin{tabular}{l|c|c|c|c}
case & $\mathrm{max}|u(p_0)|$ & $\mathrm{max}|u(p_1)|$ & $\mathrm{max}|\sigma(p_0)|$ & $\mathrm{max}|\sigma(p_1)|$ \\
\hline
"""


def gather(pressure, function_space, label):
    table = TABLE
    features = {}
    outdir = step2dir / f"step2_{label}"
    outdir.mkdir(exist_ok=True)

    for casedir in step2dir.iterdir():
        if "case" not in casedir.stem:
            continue
        case_nr = int(casedir.stem.lstrip("case"))

        for jobdir in casedir.iterdir():
            print(jobdir)
            parameter_file = jobdir / "parameters.json"
            if not parameter_file.is_file():
                continue

            params = json.loads(parameter_file.read_text())

            if params["alpha_m"] == params["alpha_m"] == 0:
                # This using newmark beta method but we are using the
                # generalized alpha method
                continue

            if params.get("function_space", "P_1") != function_space:
                continue

            if params.get("pressure") != pressure:
                continue

            if not (jobdir / "componentwise_displacement_up0.npy").is_file():
                continue

            print("Found results!")
            # loader = cardiac_benchmark.postprocess.DataLoader(jobdir / "result.h5")

            up0 = np.load(jobdir / "componentwise_displacement_up0.npy")
            up1 = np.load(jobdir / "componentwise_displacement_up1.npy")

            sp0 = np.load(jobdir / "von_Mises_stress_sp0.npy")
            sp1 = np.load(jobdir / "von_Mises_stress_sp1.npy")
            features[case_nr] = {
                "max_up0": np.linalg.norm(up0, axis=1).max(),
                "max_up1": np.linalg.norm(up1, axis=1).max(),
                "max_sp0": np.abs(sp0).max(),
                "max_sp1": np.abs(sp1).max(),
            }

            table += (
                f"{case_nr} & "
                + " & ".join(f"{v:.3f}" for v in features[case_nr].values())
                + r"\\"
                + "\n"
            )

    table += r"\end{tabular}"
    (outdir / "features.text").write_text(table)
    (outdir / "all_features.json").write_text(json.dumps(features, indent=2))
    shutil.make_archive(base_name="step2", format="zip", root_dir=outdir)
    shutil.copy("step2.zip", step2dir / f"step2_{label}.zip")


def main():

    data = [
        ("bestel", "P_1", "P1"),
        ("none", "P_1", "zero_pressure_P1"),
        ("bestel", "P_2", "P2"),
        ("none", "P_2", "zero_pressure_P2"),
    ]
    for pressure, function_space, label in data:
        gather(pressure, function_space, label)


if __name__ == "__main__":
    main()
