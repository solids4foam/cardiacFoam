#----------------------------------------------------------------------------#
# Module
#     summarize
#
# Description
#     Tables of the leadFieldECG studies from the verifier summaries
#     postProcessing/manufacturedLeadFieldECGSummary_ECG.dat of the run copies
#     in one or more results folders. Run copies are grouped by their name
#     with the token N<n> removed; the rest names the configuration. For every
#     configuration it prints, separately for the far torso electrodes
#     (names torso*) and the torso electrodes next to the heart (near*), as
#     the maximum over the group of the Linf over samples:
#       - the reciprocity gap: leadField-direct and
#         leadFieldStatic-directStatic, relative to the largest |exact|;
#       - the static errors against the exact potential at the cell
#         centroid, for the lead field and the direct solve;
#     with the order between consecutive levels of h = (cells/3)^(-1/3)
#     (1/N on hexahedra). For the layered reference it also prints:
#       - the static errors against the exact potential at the electrode
#         point;
#       - with run copies <name>_<mesh> and <name>Alpha0_<mesh>, the error
#         of the alpha part: the difference of the two static channels minus
#         the difference of their exact centroid potentials, for the left
#         (names *Left) and right (*Right) electrodes;
#       - with run copies layeredAlpha0_<mesh> and slab_<mesh>, the
#         difference of their ECGs torsoRight - torsoLeft, from the channel
#         on Vm_h and from the static channel, and each ECG's error against
#         the exact one.
#
# Usage
#     python3 summarize.py leadFieldConvergence/results
#     python3 summarize.py leadFieldConvergence/results leadFieldReciprocity/results
#     python3 summarize.py leadFieldOrganJump/results
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

SUMMARY = "postProcessing/manufacturedLeadFieldECGSummary_ECG.dat"
SAMPLES = "postProcessing/manufacturedLeadFieldECG_ECG.dat"
N_RE = re.compile(r"(?:^|_)N(\d+)(?=_|$)")

GROUPS = ("torso", "near")

STATIC = (
    ("leadFieldStatic-exactCentroid", "torso", "LF torso"),
    ("directStatic-exactCentroid", "torso", "direct torso"),
    ("leadFieldStatic-exactCentroid", "near", "LF near"),
    ("directStatic-exactCentroid", "near", "direct near"),
)

STATIC_POINT = (
    ("leadFieldStatic-exactPoint", "torso", "LF torso"),
    ("directStatic-exactPoint", "torso", "direct torso"),
    ("leadFieldStatic-exactPoint", "near", "LF near"),
    ("directStatic-exactPoint", "near", "direct near"),
)


def read_summary(path: Path) -> tuple[dict[str, float], dict[str, list[float]]]:
    header, rows = {}, {}
    for line in path.read_text().splitlines():
        words = line.split()
        if line.startswith("#"):
            if len(words) == 3:
                try:
                    header[words[1]] = float(words[2])
                except ValueError:
                    header[words[1]] = words[2]
        elif len(words) == 4:
            rows[words[0]] = [float(w) for w in words[1:]]
    return header, rows


def load(folders: list[Path]) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for folder in folders:
        for run in sorted(folder.iterdir()):
            summary = run / SUMMARY
            n = N_RE.search(run.name)
            if not summary.exists() or not n:
                continue
            header, rows = read_summary(summary)
            groups.setdefault(N_RE.sub("", run.name).strip("_"), []).append({
                "N": int(n.group(1)),
                "cells": int(header["cells"]),
                "h": (header["cells"]/3)**(-1/3),
                "scale": header["exactScale"],
                "reference": header["reference"],
                "rows": rows,
                "path": run,
            })
    for runs in groups.values():
        runs.sort(key=lambda r: r["N"])
    return groups


def order(coarse: float, fine: float, hCoarse: float, hFine: float) -> str:
    if coarse <= 0 or fine <= 0:
        return ""
    return f"{math.log(coarse/fine)/math.log(hCoarse/hFine):.2f}"


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|"*len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def group_max(rows: dict[str, list[float]], quantity: str, group: str) -> float:
    """Largest Linf of quantity over the electrodes of the group."""
    prefix = quantity + "_" + group
    values = [v[2] for k, v in rows.items() if k.startswith(prefix)]
    if not values:
        raise SystemExit(f"no rows {prefix}*")
    return max(values)


def reciprocity_table(runs: list[dict]) -> str:
    header = ["N", "cells"]
    for group in GROUPS:
        header += [f"leadField-direct {group}", f"static {group}"]
    rows = []
    for r in runs:
        cells = [str(r["N"]), str(r["cells"])]
        for group in GROUPS:
            for quantity in ("leadField-direct", "leadFieldStatic-directStatic"):
                cells.append(f"{group_max(r['rows'], quantity, group)/r['scale']:.3e}")
        rows.append(cells)
    return table(header, rows)


def static_table(runs: list[dict], quantities=STATIC) -> str:
    header = ["N", "h"]
    for _, _, label in quantities:
        header += [label, "order"]
    rows, previous = [], None
    for r in runs:
        cells = [str(r["N"]), f"{r['h']:.4g}"]
        for quantity, group, _ in quantities:
            e = group_max(r["rows"], quantity, group)
            cells.append(f"{e:.4e}")
            cells.append(
                order(group_max(previous["rows"], quantity, group), e,
                      previous["h"], r["h"])
                if previous else ""
            )
        rows.append(cells)
        previous = r
    return table(header, rows)


def read_samples(run: Path) -> list[dict[str, float]]:
    """Rows of the per-sample file, as column name -> value."""
    lines = (run / SAMPLES).read_text().splitlines()
    names = lines[0].lstrip("#").split()
    return [dict(zip(names, map(float, line.split()))) for line in lines[1:]]


def ordered_table(header: list[str], runs: list[dict],
                  values: list[list[float]]) -> str:
    """Values per run with the order between consecutive levels."""
    full = ["N", "h"]
    for label in header:
        full += [label, "order"]
    rows = []
    for i, r in enumerate(runs):
        cells = [str(r["N"]), f"{r['h']:.4g}"]
        for j, e in enumerate(values[i]):
            cells.append(f"{e:.4e}")
            cells.append(order(values[i - 1][j], e, runs[i - 1]["h"], r["h"])
                         if i else "")
        rows.append(cells)
    return table(full, rows)


def paired(a: list[dict], b: list[dict]) -> list[tuple[dict, dict]]:
    byN = {r["N"]: r for r in b}
    return [(r, byN[r["N"]]) for r in a if r["N"] in byN]


def alpha_part_table(withAlpha: list[dict], zero: list[dict]) -> str:
    """Error of the alpha part, left and right electrodes."""
    pairs = paired(withAlpha, zero)
    values = []
    for a, z in pairs:
        sa, sz = read_samples(a["path"]), read_samples(z["path"])
        row = []
        for channel in ("leadFieldStatic", "directStatic"):
            for side in ("Left", "Right"):
                e = 0.0
                for ra, rz in zip(sa, sz):
                    for name in ra:
                        if not (name.startswith(channel + "_")
                                and name.endswith(side)):
                            continue
                        electrode = name.split("_", 1)[1]
                        exact = "exactCentroid_" + electrode
                        e = max(e, abs((ra[name] - rz[name])
                                       - (ra[exact] - rz[exact])))
                row.append(e)
        values.append(row)
    return ordered_table(["LF left", "LF right", "direct left", "direct right"],
                         [a for a, _ in pairs], values)


def invariance_table(zero: list[dict], slab: list[dict]) -> str:
    """ECG torsoRight - torsoLeft of layered alpha 0 against the slab."""
    header = ["N", "LF difference", "static difference",
              "LF error layered", "LF error slab"]
    rows = []
    for z, s in paired(zero, slab):
        sz, ss = read_samples(z["path"]), read_samples(s["path"])

        def ecg(row: dict, channel: str) -> float:
            return row[channel + "_torsoRight"] - row[channel + "_torsoLeft"]

        dLF = max(abs(ecg(a, "leadField") - ecg(b, "leadField"))
                  for a, b in zip(sz, ss))
        dStatic = max(abs(ecg(a, "leadFieldStatic") - ecg(b, "leadFieldStatic"))
                      for a, b in zip(sz, ss))
        eZ = max(abs(ecg(a, "leadField") - ecg(a, "exactCentroid")) for a in sz)
        eS = max(abs(ecg(b, "leadField") - ecg(b, "exactCentroid")) for b in ss)
        rows.append([str(z["N"]), f"{dLF:.3e}", f"{dStatic:.3e}",
                     f"{eZ:.4e}", f"{eS:.4e}"])
    return table(header, rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="+",
                        help="folders of run copies")
    args = parser.parse_args()
    groups = load(args.results)
    if not groups:
        raise SystemExit("no run copies with a lead-field verifier summary found")
    for config in sorted(groups):
        print(f"\n### {config}\n")
        print("reciprocity gap relative to the largest |exact|:\n")
        print(reciprocity_table(groups[config]))
        print("\nstatic errors against the exact potential:\n")
        print(static_table(groups[config]))
        if groups[config][0]["reference"] == "layered":
            print("\nstatic errors against the exact potential at the "
                  "electrode point:\n")
            print(static_table(groups[config], STATIC_POINT))
        problem, _, mesh = config.partition("_")
        zero = f"{problem}Alpha0_{mesh}"
        if zero in groups:
            print(f"\nerror of the alpha part (against {zero}):\n")
            print(alpha_part_table(groups[config], groups[zero]))
        if problem == "layeredAlpha0" and f"slab_{mesh}" in groups:
            print("\nECG torsoRight - torsoLeft against slab:\n")
            print(invariance_table(groups[config], groups[f"slab_{mesh}"]))


if __name__ == "__main__":
    main()
