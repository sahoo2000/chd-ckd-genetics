#!/usr/bin/env python3
"""
07_summarize.py — parse regenie output, check calibration, flag unstable estimates.

Built around the failure modes that produced spurious results before:
  * regenie reports LOG10P, not P. Converting wrongly is the classic error.
  * singleton masks with <20 carriers give |BETA| > 5 and meaningless ORs.
  * lambda_GC is deflated for sparse burden statistics — that is expected, not a bug.
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd
from scipy.stats import chi2

MIN_CARRIERS = 20   # below this, a burden estimate is not interpretable


def load(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep=r"\s+", comment="#", engine="python")
        if "LOG10P" not in df.columns:
            print(f"[skip] {p}: no LOG10P column")
            continue
        df["source"] = os.path.basename(p)
        frames.append(df)
    if not frames:
        sys.exit("ERROR: no parsable regenie output found.")
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="burden_*.regenie")
    ap.add_argument("--n-total", type=float, default=450000,
                    help="analysis sample size, for carrier-count estimates")
    ap.add_argument("--out", default="summary.tsv")
    a = ap.parse_args()

    files = sorted(glob.glob(a.glob))
    if not files:
        sys.exit(f"ERROR: nothing matched {a.glob}")
    print(f"[read] {len(files)} file(s)")
    d = load(files)

    d["P"] = 10.0 ** (-d["LOG10P"])
    if "ID" in d.columns:
        parts = d["ID"].astype(str).str.split(".", n=1, expand=True)
        d["SET"] = parts[0]
        d["MASK"] = parts[1] if parts.shape[1] > 1 else ""
    n = a.n_total
    d["carriers_est"] = (d.get("A1FREQ", np.nan).astype(float) * 2 * n).round(0)
    d["unstable"] = (d.carriers_est < MIN_CARRIERS) | (d.get("BETA", pd.Series(np.nan)).abs() > 5)

    ntests = len(d)
    nsets = d["SET"].nunique() if "SET" in d else ntests
    bonf_all = 0.05 / max(ntests, 1)
    bonf_set = 0.05 / max(nsets, 1)

    print("\n" + "=" * 68)
    print(f"tests: {ntests:,}   distinct sets: {nsets:,}")
    print(f"Bonferroni  all tests: {bonf_all:.3g}   per set: {bonf_set:.3g}")
    if nsets == 1:
        print("single pooled set -> alpha = 0.05 is the correct threshold")
    print("=" * 68)

    # ---- calibration ----
    print("\nCALIBRATION")
    for thr in (0.05, 0.01, 1e-3):
        obs = int((d.P < thr).sum()); exp = thr * ntests
        print(f"  p<{thr:<7g} observed {obs:>6,}  expected {exp:>8.1f}  ratio {obs/max(exp,1e-9):.2f}")
    valid = d.P.dropna()
    valid = valid[(valid > 0) & (valid < 1)]
    if len(valid) > 10:
        lam = np.median(chi2.isf(valid, 1)) / chi2.ppf(0.5, 1)
        print(f"  lambda_GC {lam:.3f}"
              + ("  (deflation is expected for sparse burden tests)" if lam < 0.9 else ""))

    # ---- unstable ----
    nu = int(d.unstable.sum())
    if nu:
        print(f"\nUNSTABLE ESTIMATES: {nu} test(s) with <{MIN_CARRIERS} carriers or |BETA|>5")
        print("  These are sparse-data separation artefacts. Do not report them.")
        cols = [c for c in ("SET", "MASK", "TEST", "A1FREQ", "carriers_est", "BETA", "SE", "P")
                if c in d.columns]
        print(d[d.unstable].nsmallest(min(nu, 8), "P")[cols]
              .to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # ---- results ----
    good = d[~d.unstable].copy()
    good = good.sort_values("P")
    if len(good):
        # Benjamini-Hochberg: q_(i) = min over j>=i of (n/j * p_(j)) -> cummin runs
        # from the LARGEST p-value backwards, not forwards.
        raw = good.P.values * len(good) / (np.arange(len(good)) + 1)
        good["q_BH"] = np.minimum.accumulate(raw[::-1])[::-1].clip(max=1)
        thr = 0.05 if nsets == 1 else bonf_set
        sig = good[good.P < thr]
        print(f"\nRESULTS (stable estimates only, n={len(good)})")
        print(f"  passing alpha={thr:.3g}: {len(sig)}")
        cols = [c for c in ("SET", "MASK", "TEST", "A1FREQ", "carriers_est",
                            "BETA", "SE", "P", "q_BH", "source") if c in good.columns]
        print(good.head(15)[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
        good.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
        print(f"\n[write] {a.out}")
    else:
        print("\nNo stable estimates survived filtering.")


if __name__ == "__main__":
    main()
