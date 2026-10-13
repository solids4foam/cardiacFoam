#!/usr/bin/env python3
"""Independent surface-integral and common-sampling checks for eikonal ECG."""
import argparse
import json
import math
from pathlib import Path

from checkDualECG import cohort, series


def gauss_rule(n):
    """Gauss-Legendre rule on [0,1], independent of the C++ reference code."""
    result = []
    for i in range(n):
        x = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(64):
            previous, value = 1.0, x
            for j in range(2, n + 1):
                previous, value = value, ((2*j-1)*x*value-(j-1)*previous)/j
            derivative = n*(x*value-previous)/(x*x-1)
            delta = value/derivative
            x -= delta
            if abs(delta) < 1e-14:
                break
        else:
            raise AssertionError("Legendre root did not converge")
        result.append(((x+1)/2, 1/((1-x*x)*derivative*derivative)))
    return result


def boundary_coefficients(n):
    # Constant L=(1,2,3): sum of face voltage differences, with identical
    # faces by permutation symmetry. No numerical voltage gradient is used.
    rule = gauss_rule(n)
    cosine = sine = 0.0
    for x, wx in rule:
        for y, wy in rule:
            low = 2*math.pi*math.exp(x+y)
            high = 2*math.pi*math.exp(1+x+y)
            cosine += 6*wx*wy*(math.sin(low)-math.sin(high))
            sine += 6*wx*wy*(math.cos(high)-math.cos(low))
    return cosine, sine


def check(case):
    out = case / "postProcessing"
    lh, lr = series(out / "manufacturedEikonalLeadFieldECG_LeadField.dat")
    yh, yr = series(out / "eikonalLeadFieldECG.dat")
    ph, pr = series(out / "manufacturedEikonalECG_ECG.dat")
    zh, zr = series(out / "eikonalECG.dat")
    assert len(lr) == len(yr) == len(pr) == len(zr) == 201
    for records in zip(lr, yr, pr, zr):
        assert max(r[0] for r in records)-min(r[0] for r in records) < 1e-12
    assert lr[0][0] == 0 and lr[-1][0] == 2
    for lead in ("positive", "negative"):
        assert all(a[lh.index(lead+"_numeric")] == b[yh.index(lead)]
                   for a, b in zip(lr, yr))
    for electrode in ("E1", "E2", "E3", "E4", "E5"):
        assert all(a[ph.index("numeric_"+electrode)] == b[zh.index(electrode)]
                   for a, b in zip(pr, zr))
    c, s = boundary_coefficients(80)
    reference_mismatch = quadrature_delta = 0.0
    for row in lr:
        exact = c*math.cos(2*math.pi*row[0])+s*math.sin(2*math.pi*row[0])
        for lead, sign in (("positive", 1), ("negative", -1)):
            reference_mismatch = max(reference_mismatch,
                abs(row[lh.index(lead+"_exact")]-sign*exact))
            quadrature_delta = max(quadrature_delta,
                row[lh.index(lead+"_quadratureDelta")])
        assert abs(row[lh.index("positive_numeric")]
                   +row[lh.index("negative_numeric")]) < 1e-12
    assert reference_mismatch < 5e-9, "Independent surface reference mismatch"
    leads = cohort(lh, lr, [(lead+"_numeric", lead+"_exact")
                           for lead in ("positive", "negative")])
    pseudo = cohort(ph, pr, [("numeric_"+e, "refQ48_"+e)
                            for e in ("E3", "E4", "E5")])
    assert leads["relative_RMS_error"] < 0.05, "Reciprocal functional error bound"
    assert pseudo["relative_RMS_error"] < 0.1, "Pseudo functional error bound"
    assert quadrature_delta < max(1e-13, 1e-3*leads["RMS_error"])
    source = {}
    for line in (out / "manufacturedEikonalActivationTime.dat").read_text().splitlines():
        columns = line.split()
        if columns and columns[0] in ("activationTime", "gradActivationTime", "gradReconActivationTime"):
            source[columns[0]] = list(map(float, columns[1:]))
    assert all(math.isfinite(v) for values in source.values() for v in values)
    assert source["activationTime"][1] < 3e-4 and source["activationTime"][2] < 2e-3
    assert source["gradActivationTime"][1] < 1e-2
    assert source["gradReconActivationTime"][1] < 2e-3
    fixture_log = (case / "log.setManufacturedLeadFields").read_text()
    for template in ("manufactured", "endo", "mid", "epi"):
        assert "PASS: affine "+template+" template" in fixture_log
    result = {"samples": len(lr), "leadFieldECG": leads, "pseudoECG_raw": pseudo,
              "surface64_96_max_delta": quadrature_delta,
              "independent_surface80_max_mismatch": reference_mismatch,
              "source_L1_L2_Linf": source,
              "comparison": "Separate sensitivities and independent references; no model-equivalence claim."}
    (out / "eikonalDualECGSummary.json").write_text(json.dumps(result, indent=2)+"\n")
    print("PASS: eikonal dual ECG, four template controls and independent surface reference",
          json.dumps(result))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    check(parser.parse_args().case)
