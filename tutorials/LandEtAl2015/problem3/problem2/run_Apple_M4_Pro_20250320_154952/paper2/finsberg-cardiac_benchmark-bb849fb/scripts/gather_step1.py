import shutil
from pathlib import Path

import numpy as np

import cardiac_benchmark

here = Path(__file__).absolute().parent


def main(jobid, label):
    resultdir = Path(f"/global/D1/homes/henriknf/cardiac_benchmark/step1/{jobid}")
    disp_path = here / ".." / "data" / "displacement_points.npz"
    vol_path = here / ".." / "data" / "computed_vols.npz"

    loader = cardiac_benchmark.postprocess.DataLoader(resultdir / "result.h5")

    time = loader.time_stamps
    up0 = np.load(resultdir / "componentwise_displacement_up0.npy")
    up1 = np.load(resultdir / "componentwise_displacement_up1.npy")
    vol = np.load(resultdir / "volume.npy")

    vol_true = cardiac_benchmark.postprocess.load_true_volume_data(vol_path=vol_path)
    disp_true = cardiac_benchmark.postprocess.load_true_displacement_data(
        disp_path=disp_path,
    )

    outdir = Path(f"step1_{label}")
    outdir.mkdir(exist_ok=True)

    cardiac_benchmark.postprocess.plot_componentwise_displacement_comparison(
        up0,
        up1,
        time,
        disp_true.up0,
        disp_true.up1,
        disp_true.time,
        fname=outdir / "componentwise_displacement_comparison.png",
    )

    cardiac_benchmark.postprocess.plot_volume_comparison(
        vol,
        time,
        vol_true,
        disp_true.time,
        fname=outdir / "volume_comparison.png",
    )

    np.save(outdir / "up0.npy", up0)
    np.save(outdir / "up1.npy", up1)
    np.save(outdir / "time.npy", time)
    np.save(outdir / "vol.npy", vol)
    shutil.make_archive(base_name="step1", format="zip", root_dir=outdir)
    shutil.copy("step1.zip", resultdir.parent / f"step1_{label}.zip")


if __name__ == "__main__":
    data = {
        491494: "P1",
        491461: "zero_pressure_P1",
        491630: "P2",
        491631: "zero_pressure_P2",
    }
    for jobid, label in data.items():
        main(jobid, label)
