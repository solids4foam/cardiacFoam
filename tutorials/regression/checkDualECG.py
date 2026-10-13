#!/usr/bin/env python3
"""Independent common-cohort error/finite-value/quadrature checks for dual ECG MMS."""
import argparse,json,math
from pathlib import Path

def series(path):
    lines=path.read_text().splitlines()
    header=next(l[1:].split() for l in lines if l.startswith("#") and "time" in l.split())
    rows=[list(map(float,l.split())) for l in lines if l.strip() and not l.startswith("#")]
    assert rows and all(len(r)==len(header) for r in rows),(path,header)
    assert all(math.isfinite(v) for r in rows for v in r),path
    assert all(rows[i][0]<rows[i+1][0] for i in range(len(rows)-1)),path
    return header,rows

def cohort(header,rows,pairs):
    error=reference=0.
    for row in rows:
        for numeric,exact in pairs:
            y=row[header.index(numeric)]; e=row[header.index(exact)]
            error+=(y-e)**2; reference+=e**2
    assert reference>1e-20,"Analytically zero cohort cannot define relative error"
    count=len(rows)*len(pairs)
    return {"relative_RMS_error":math.sqrt(error/reference),"RMS_error":math.sqrt(error/count),"RMS_reference":math.sqrt(reference/count)}

def check(case):
    out=case/"postProcessing"
    lh,lr=series(out/"manufacturedLeadFieldECG_ECG.dat")
    ph,pr=series(out/"manufacturedPseudoECG_Pseudo.dat")
    assert len(lr)==len(pr),"Observer reference sample counts differ"
    assert all(abs(a[0]-b[0])<1e-12 for a,b in zip(lr,pr)),"Observer reference times differ"
    leads=cohort(lh,lr,[("positive_numeric","positive_exact"),("negative_numeric","negative_exact")])
    pseudo=cohort(ph,pr,[("numeric_"+e,"refQ96_"+e) for e in ["E3","E4","E5"]])
    delta=sum(r[ph.index("deltaQuadratureQ48_Q96_"+e)]**2 for r in pr for e in ["E3","E4","E5"])
    delta=math.sqrt(delta/(3*len(pr)))
    assert delta<=max(1e-14,1e-3*pseudo["RMS_error"]),"Quadrature uncertainty is not subordinate"
    assert delta/pseudo["RMS_reference"]<1e-8,"Quadrature reference control failed"
    assert leads["relative_RMS_error"]<.1 and pseudo["relative_RMS_error"]<.4,"Functional error bound failed"
    result={"samples":len(lr),"leadFieldECG":leads,"pseudoECG_raw":pseudo,"pseudo_per_electrode":{e:cohort(ph,pr,[("numeric_"+e,"refQ96_"+e)]) for e in ["E3","E4","E5"]},"quadrature48_96_RMS_delta":delta,"quadrature_relative_RMS_delta":delta/pseudo["RMS_reference"],"comparison":"Different functionals, compared only through independently normalized error cohorts; pseudo sigmaExtracellular 0 unnormalised case functional"}
    (out/"dualECGSummary.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: dual observers and subordinate quadrature",json.dumps(result))
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("case",type=Path);args=parser.parse_args();check(args.case)
