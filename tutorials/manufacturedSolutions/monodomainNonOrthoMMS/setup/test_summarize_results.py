from pathlib import Path

import pytest

from summarize_results import (
    build_arg_parser,
    build_summary_rows,
    convergence_order,
    parse_dat_file,
)


def test_convergence_order_second_order_data():
    # errors halving-squared each refinement -> order 2 exactly
    errors = [8.0e-2, 2.0e-2, 5.0e-3, 1.25e-3]
    orders = convergence_order(errors)
    assert orders[0] is None
    for order in orders[1:]:
        assert order == pytest.approx(2.0, abs=1e-9)


def test_parse_dat_file(tmp_path):
    dat = tmp_path / "3D_10_cells_implicit.dat"
    dat.write_text(
        "Manufactured-solution error summary (t = 0.2):\n"
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     1.0e-02   2.0e-02   3.0e-02\n"
        "u1     4.0e-02   5.0e-02   6.0e-02\n"
        "-------------------------------------------------\n\n"
        "Simulation summary:\n"
        "-------------------\n"
        "Number of cells (N)   = 10\n"
        "Grid spacing (dx)     = 0.1\n"
        "Time step (dt)        = 0.00892857\n"
        "-------------------\n"
    )
    summary = parse_dat_file(dat)
    assert summary.dx == 0.1
    assert summary.dt == 0.00892857
    assert summary.l2_vm == 2.0e-02


def test_build_summary_rows_end_to_end(tmp_path):
    results_dir = tmp_path
    errors = {10: 8.0e-2, 20: 2.0e-2, 40: 5.0e-3, 80: 1.25e-3}
    for n, err in errors.items():
        case_dir = results_dir / "0.0"
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / f"3D_{n}_cells_implicit.dat").write_text(
            "Manufactured-solution error summary (t = 0.2):\n"
            "Field     L1-error       L2-error       Linf-error\n"
            f"Vm     {err}   {err}   {err}\n"
            "-------------------------------------------------\n\n"
            "Simulation summary:\n"
            f"Grid spacing (dx)     = {1.0 / n}\n"
            "Time step (dt)        = 0.001\n"
        )
        (case_dir / f"log.checkMesh.{n}").write_text(
            "Mesh non-orthogonality Max: 0.01 average: 0.001\nMesh OK.\n"
        )

    rows = build_summary_rows(results_dir, amplitudes=["0.0"], resolutions=[10, 20, 40, 80])
    assert len(rows) == 4
    assert rows[0]["order_Vm"] is None
    assert rows[-1]["order_Vm"] == 2.0


def test_cli_amplitudes_preserve_trailing_zero():
    # Regression test for the original bug: run_nonortho_sweep.sh names
    # result directories with literal shell tokens like "0.10". The CLI
    # used to declare --amplitudes with type=float, so argparse itself
    # mangled "0.10" into the float 0.1 (str(float("0.10")) == "0.1")
    # before build_summary_rows ever saw it, causing a directory-lookup
    # mismatch and a crash. This test exercises the actual argparse layer
    # where that bug lived: reverting --amplitudes back to type=float
    # makes this assertion fail immediately (args.amplitudes would be
    # [0.1, 0.2], not ["0.10", "0.20"]).
    parser = build_arg_parser()
    args = parser.parse_args(
        ["some_dir", "--amplitudes", "0.10", "0.20", "--resolutions", "10"]
    )
    assert args.amplitudes == ["0.10", "0.20"]


def test_build_summary_rows_accepts_trailing_zero_amplitude_label(tmp_path):
    # Not a regression test for the argparse bug above (build_summary_rows
    # never mangles its inputs -- it just does path joins on whatever
    # string it's given, whether or not that string is "float-safe").
    # This checks the separate, still-legitimate concern that a
    # "0.10"-named results directory is looked up and parsed without
    # error.
    results_dir = tmp_path
    case_dir = results_dir / "0.10"
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "3D_10_cells_implicit.dat").write_text(
        "Manufactured-solution error summary (t = 0.2):\n"
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     8.0e-02   8.0e-02   8.0e-02\n"
        "-------------------------------------------------\n\n"
        "Simulation summary:\n"
        "Grid spacing (dx)     = 0.1\n"
        "Time step (dt)        = 0.00892857\n"
    )
    (case_dir / "log.checkMesh.10").write_text(
        "Mesh non-orthogonality Max: 0.01 average: 0.001\nMesh OK.\n"
    )

    rows = build_summary_rows(results_dir, amplitudes=["0.10"], resolutions=[10])
    assert len(rows) == 1
    assert rows[0]["A"] == "0.10"
    assert rows[0]["order_Vm"] is None
