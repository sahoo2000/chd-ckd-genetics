#!/usr/bin/env python3
"""
01_derive_phenotypes.py — kidney-function traits + congenital case definitions from UK Biobank.

Produces the phenotype file regenie needs, with eGFR as the PRIMARY outcome. The whole
point of the design is to stop conditioning on a rare binary comorbidity (~27 people)
and use a continuous kidney trait measured in essentially everyone (~450,000).

Inputs
  --tab       UKB main dataset, tab-delimited, with the field columns named f.<field>.<inst>.<arr>
              (the default naming from ukbconv / ukbtools). Adjust FIELD_COL if yours differ.
  --hesin     optional: HES diagnoses long file (eid, diag_icd10) for ICD ascertainment.
              If omitted, ICD codes are read from the first-occurrence array field 41270.

Outputs
  pheno_quant.tsv   FID IID eGFR_cr eGFR_cys eGFR_cr_cys  (RINT-transformed columns too)
  pheno_binary.tsv  FID IID CHD CAKUT CKD CHD_CAKUT
  pheno_summary.txt counts and distributions — READ THIS before running regenie

Equations: CKD-EPI 2021 race-free (Inker et al., NEJM 2021;385:1737). Verify the constants
against the paper before publishing — they are transcribed here, not derived.
"""
import argparse, sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- UKB field IDs
F_CREAT   = 30700   # Creatinine, serum (umol/L)
F_CYS     = 30720   # Cystatin C (mg/L)
F_SEX     = 31      # Sex (0=Female, 1=Male)
F_AGE     = 21022   # Age at recruitment
F_ICD10   = 41270   # Diagnoses - ICD10 (array)

CREAT_UMOL_PER_MGDL = 88.4

# ICD-10 prefixes
CHD_CODES   = tuple(f"Q{n}" for n in range(20, 27))   # Q20-Q26 congenital heart / great vessels
CAKUT_CODES = tuple(f"Q{n}" for n in range(60, 65))   # Q60-Q64 kidney & urinary tract
CKD_CODES   = ("N18",)                                 # chronic kidney disease


def field_cols(df, field, inst=0):
    """Return columns for a UKB field at a given instance, across all array indices."""
    pre = f"f.{field}.{inst}."
    cols = [c for c in df.columns if c.startswith(pre)]
    if not cols:                       # tolerate the 30700-0.0 style too
        pre2 = f"{field}-{inst}."
        cols = [c for c in df.columns if c.startswith(pre2)]
    return cols


def first_col(df, field, inst=0):
    c = field_cols(df, field, inst)
    if not c:
        sys.exit(f"ERROR: field {field} (instance {inst}) not found in the input table.")
    return df[c[0]]


def egfr_cr(scr_mgdl, age, female):
    """CKD-EPI 2021 creatinine, race-free."""
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.241, -0.302)
    r = scr_mgdl / kappa
    return (142.0
            * np.minimum(r, 1.0) ** alpha
            * np.maximum(r, 1.0) ** -1.200
            * 0.9938 ** age
            * np.where(female, 1.012, 1.0))


def egfr_cys(scys, age, female):
    """CKD-EPI cystatin C (2012), race-free."""
    r = scys / 0.8
    return (133.0
            * np.minimum(r, 1.0) ** -0.499
            * np.maximum(r, 1.0) ** -1.328
            * 0.996 ** age
            * np.where(female, 0.932, 1.0))


def egfr_cr_cys(scr_mgdl, scys, age, female):
    """CKD-EPI 2021 combined creatinine-cystatin C, race-free."""
    kappa = np.where(female, 0.7, 0.9)
    alpha = np.where(female, -0.219, -0.144)
    rc = scr_mgdl / kappa
    ry = scys / 0.8
    return (135.0
            * np.minimum(rc, 1.0) ** alpha
            * np.maximum(rc, 1.0) ** -0.544
            * np.minimum(ry, 1.0) ** -0.323
            * np.maximum(ry, 1.0) ** -0.778
            * 0.9961 ** age
            * np.where(female, 0.963, 1.0))


def rint(x):
    """Rank-based inverse normal transform, NaN-safe."""
    from scipy.stats import norm
    x = pd.Series(x)
    ok = x.notna()
    r = x[ok].rank(method="average")
    out = pd.Series(np.nan, index=x.index)
    out[ok] = norm.ppf((r - 0.5) / ok.sum())
    return out.values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tab", required=True, help="UKB main dataset (tab-delimited)")
    ap.add_argument("--hesin", help="optional HES diagnoses long file (eid, diag_icd10)")
    ap.add_argument("--out-prefix", default="pheno")
    ap.add_argument("--winsor-sd", type=float, default=5.0,
                    help="winsorise eGFR beyond this many SD (0 disables)")
    a = ap.parse_args()

    print(f"[read] {a.tab}", flush=True)
    df = pd.read_csv(a.tab, sep="\t", low_memory=False)
    eid = df.columns[0]
    print(f"[read] {len(df):,} participants")

    sex_male = first_col(df, F_SEX).astype(float)          # 1 = male
    female = (sex_male == 0).values
    age = first_col(df, F_AGE).astype(float).values

    scr_umol = first_col(df, F_CREAT).astype(float).values
    scr = scr_umol / CREAT_UMOL_PER_MGDL
    scys = first_col(df, F_CYS).astype(float).values

    q = pd.DataFrame({"FID": df[eid], "IID": df[eid]})
    q["eGFR_cr"] = egfr_cr(scr, age, female)
    q["eGFR_cys"] = egfr_cys(scys, age, female)
    q["eGFR_cr_cys"] = egfr_cr_cys(scr, scys, age, female)

    # winsorise then RINT — regenie's --apply-rint does this too, but doing it here
    # makes the distribution auditable before it reaches the model
    for c in ["eGFR_cr", "eGFR_cys", "eGFR_cr_cys"]:
        v = q[c].values.astype(float)
        if a.winsor_sd > 0:
            m, s = np.nanmean(v), np.nanstd(v)
            v = np.clip(v, m - a.winsor_sd * s, m + a.winsor_sd * s)
            q[c] = v
        q[c + "_rint"] = rint(v)

    # ---------------------------------------------------------------- ICD codes
    print("[icd] ascertaining congenital diagnoses", flush=True)
    if a.hesin:
        h = pd.read_csv(a.hesin, sep="\t", low_memory=False,
                        usecols=lambda c: c.lower() in ("eid", "diag_icd10"))
        h.columns = [c.lower() for c in h.columns]
        codes = h.groupby("eid")["diag_icd10"].apply(lambda s: set(s.dropna().astype(str)))
    else:
        icd_cols = field_cols(df, F_ICD10)
        if not icd_cols:
            sys.exit("ERROR: no ICD source. Provide --hesin or include field 41270.")
        sub = df[[eid] + icd_cols].set_index(eid)
        codes = sub.apply(lambda r: set(str(v) for v in r.dropna()), axis=1)

    def flag(prefixes):
        s = codes.apply(lambda cs: int(any(c.startswith(prefixes) for c in cs)))
        return s.reindex(df[eid]).fillna(0).astype(int).values

    b = pd.DataFrame({"FID": df[eid], "IID": df[eid]})
    b["CHD"] = flag(CHD_CODES)
    b["CAKUT"] = flag(CAKUT_CODES)
    b["CKD"] = flag(CKD_CODES)
    b["CHD_CAKUT"] = ((b.CHD == 1) & (b.CAKUT == 1)).astype(int)

    q.to_csv(f"{a.out_prefix}_quant.tsv", sep="\t", index=False, na_rep="NA")
    b.to_csv(f"{a.out_prefix}_binary.tsv", sep="\t", index=False, na_rep="NA")

    # ---------------------------------------------------------------- summary
    lines = ["UK Biobank phenotype summary", "=" * 62, f"participants: {len(df):,}", ""]
    lines.append("QUANTITATIVE (primary outcomes)")
    for c in ["eGFR_cr", "eGFR_cys", "eGFR_cr_cys"]:
        v = q[c].dropna()
        lines.append(f"  {c:<14} n={len(v):>7,}  mean={v.mean():6.1f}  sd={v.std():5.1f}  "
                     f"median={v.median():6.1f}  <60={100*(v<60).mean():4.1f}%")
    lines += ["", "BINARY (ascertainment check — NOT the primary design)"]
    for c in ["CHD", "CAKUT", "CKD", "CHD_CAKUT"]:
        n = int(b[c].sum())
        lines.append(f"  {c:<14} cases={n:>7,}  ({100*n/len(b):.3f}%)")
    n_comorbid = int(b.CHD_CAKUT.sum())
    lines += ["", "-" * 62]
    if n_comorbid < 200:
        lines.append(f"WARNING: CHD_CAKUT has {n_comorbid} cases. A burden test on this "
                     f"phenotype\n         has no power at any plausible effect size. Use the "
                     f"quantitative\n         eGFR outcomes as the primary analysis.")
    else:
        lines.append(f"CHD_CAKUT n={n_comorbid} — larger than expected; recheck the code lists.")
    txt = "\n".join(lines)
    open(f"{a.out_prefix}_summary.txt", "w").write(txt + "\n")
    print("\n" + txt)


if __name__ == "__main__":
    main()
