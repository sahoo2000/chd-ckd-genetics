#!/usr/bin/env python3
"""
06_test_prs.py — is CKD common-variant risk elevated in congenital heart disease cases?

Logistic regression of case status on the standardised CKD polygenic score, adjusted for
age, sex and PCs. Reports OR per SD with a 95% CI.

Direction note: a CKD-lowering-eGFR score raises risk, so sign the score so that HIGHER
means MORE CKD risk before interpreting. The --flip-sign flag does this when your
summary statistics are on the eGFR (higher = better kidney function) scale.
"""
import argparse, glob, sys
import numpy as np
import pandas as pd
import statsmodels.api as sm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prs-glob", required=True, help="glob for PLINK2 .sscore files")
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--covar", required=True)
    ap.add_argument("--outcomes", nargs="+", default=["CHD"])
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--flip-sign", action="store_true",
                    help="negate the score (use when weights are on the eGFR scale, so "
                         "that higher score = higher CKD risk)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ph = pd.read_csv(a.pheno, sep="\t")
    cv = pd.read_csv(a.covar, sep="\t")
    cov_cols = ["age", "age2", "sex", "age_sex"] + [f"PC{i+1}" for i in range(a.n_pcs)]
    missing = [c for c in cov_cols if c not in cv.columns]
    if missing:
        sys.exit(f"ERROR: covariate columns missing: {missing}")

    files = sorted(glob.glob(a.prs_glob))
    if not files:
        sys.exit(f"ERROR: no files matched {a.prs_glob}")

    rows = []
    for f in files:
        tag = f.split("prs_")[-1].replace(".sscore", "")
        s = pd.read_csv(f, sep="\t")
        idc = "IID" if "IID" in s.columns else s.columns[0]
        scol = next((c for c in s.columns if "SCORE" in c.upper() and "AVG" in c.upper()),
                    next((c for c in s.columns if "SCORE" in c.upper()), None))
        if scol is None:
            print(f"[skip] no score column in {f}")
            continue
        s = s[[idc, scol]].rename(columns={idc: "IID", scol: "score"})

        d = ph.merge(s, on="IID").merge(cv, on="IID")
        d = d.dropna(subset=["score"] + cov_cols)
        if a.flip_sign:
            d["score"] = -d["score"]
        d["score_z"] = (d["score"] - d["score"].mean()) / d["score"].std()

        for outcome in a.outcomes:
            if outcome not in d.columns:
                continue
            y = d[outcome].astype(float)
            ncase = int(y.sum())
            if ncase < 50:
                rows.append(dict(threshold=tag, outcome=outcome, n_case=ncase,
                                 n_total=len(d), OR=np.nan, lo=np.nan, hi=np.nan,
                                 p=np.nan, note="too few cases to fit"))
                continue
            X = sm.add_constant(d[["score_z"] + cov_cols].astype(float))
            try:
                m = sm.Logit(y, X).fit(disp=0)
                b, se, p = m.params["score_z"], m.bse["score_z"], m.pvalues["score_z"]
                rows.append(dict(threshold=tag, outcome=outcome, n_case=ncase,
                                 n_total=len(d), OR=np.exp(b),
                                 lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                                 p=p, note=""))
            except Exception as e:
                rows.append(dict(threshold=tag, outcome=outcome, n_case=ncase,
                                 n_total=len(d), OR=np.nan, lo=np.nan, hi=np.nan,
                                 p=np.nan, note=f"fit failed: {e}"))

    r = pd.DataFrame(rows)
    r.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    print(r.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\n[write] {a.out}")
    print("\nMultiple thresholds were tested — correct for that, or pre-register one.")


if __name__ == "__main__":
    main()
