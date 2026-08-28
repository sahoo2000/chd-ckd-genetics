#!/usr/bin/env python3
"""
02_make_covariates.py — covariate file and analysis-sample keep-list for regenie.

Covariates: age, age^2, sex, age*sex, assessment centre, genotyping batch, PC1-20.
Exclusions: sex-chromosome aneuploidy, heterozygosity/missingness outliers,
            sex mismatch, withdrawal list, and (optionally) relateds.

Outputs
  covariates.tsv   FID IID age age2 sex age_sex centre batch PC1..PC20
  keep.tsv         FID IID   — the analysis sample, pass to regenie --keep
  exclusions.txt   how many participants each filter removed
"""
import argparse, sys
import numpy as np
import pandas as pd

F_SEX_SELF = 31      # Sex (self-reported / registry)
F_SEX_GEN  = 22001   # Genetic sex
F_AGE      = 21022   # Age at recruitment
F_CENTRE   = 54      # Assessment centre
F_BATCH    = 22000   # Genotype measurement batch
F_PC       = 22009   # Genetic principal components (array 1..40)
F_ANEU     = 22019   # Sex chromosome aneuploidy
F_HETMISS  = 22027   # Outliers for heterozygosity or missing rate
F_INWHITE  = 22006   # Genetic ethnic grouping (Caucasian == 1)


def field_cols(df, field, inst=0):
    pre = f"f.{field}.{inst}."
    cols = [c for c in df.columns if c.startswith(pre)]
    if not cols:
        cols = [c for c in df.columns if c.startswith(f"{field}-{inst}.")]
    return cols


def first_col(df, field, inst=0, required=True):
    c = field_cols(df, field, inst)
    if not c:
        if required:
            sys.exit(f"ERROR: field {field} not found.")
        return pd.Series(np.nan, index=df.index)
    return df[c[0]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tab", required=True)
    ap.add_argument("--withdrawals", help="one eid per line")
    ap.add_argument("--related", help="one eid per line to drop (optional; regenie's "
                                      "whole-genome model handles relatedness, so "
                                      "dropping relateds is usually unnecessary)")
    ap.add_argument("--ancestry", choices=["all", "eur"], default="eur",
                    help="'eur' keeps field 22006 == 1. Run ancestries separately, "
                         "never pooled with PCs alone.")
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--out-prefix", default="")
    a = ap.parse_args()

    df = pd.read_csv(a.tab, sep="\t", low_memory=False)
    eid = df.columns[0]
    n0 = len(df)
    print(f"[read] {n0:,} participants")

    cov = pd.DataFrame({"FID": df[eid], "IID": df[eid]})
    age = first_col(df, F_AGE).astype(float)
    sex = first_col(df, F_SEX_GEN, required=False)
    if sex.isna().all():
        sex = first_col(df, F_SEX_SELF)
    sex = sex.astype(float)

    cov["age"] = age
    cov["age2"] = age ** 2
    cov["sex"] = sex
    cov["age_sex"] = age * sex
    cov["centre"] = first_col(df, F_CENTRE).astype("Int64")
    cov["batch"] = first_col(df, F_BATCH, required=False).astype("Int64")

    pcs = field_cols(df, F_PC)
    if len(pcs) < a.n_pcs:
        sys.exit(f"ERROR: found {len(pcs)} PC columns, need {a.n_pcs}.")
    for i in range(a.n_pcs):
        cov[f"PC{i+1}"] = df[pcs[i]].astype(float)

    # ---------------------------------------------------------------- exclusions
    keep = pd.Series(True, index=df.index)
    log = []

    def drop(mask, label):
        nonlocal keep
        n_before = int(keep.sum())
        keep &= ~mask.fillna(False)
        log.append(f"  {label:<44} -{n_before - int(keep.sum()):>7,}")

    aneu = first_col(df, F_ANEU, required=False)
    drop(aneu.notna() & (aneu.astype(float) == 1), "sex chromosome aneuploidy")

    hm = first_col(df, F_HETMISS, required=False)
    drop(hm.notna() & (hm.astype(float) == 1), "heterozygosity / missingness outlier")

    self_sex = first_col(df, F_SEX_SELF).astype(float)
    gen_sex = first_col(df, F_SEX_GEN, required=False).astype(float)
    drop(gen_sex.notna() & (self_sex != gen_sex), "reported vs genetic sex mismatch")

    if a.ancestry == "eur":
        w = first_col(df, F_INWHITE, required=False).astype(float)
        drop(~(w == 1), "not in genetic Caucasian grouping (22006)")

    if a.withdrawals:
        wd = set(pd.read_csv(a.withdrawals, header=None)[0].astype(str))
        drop(df[eid].astype(str).isin(wd), "consent withdrawn")

    if a.related:
        rel = set(pd.read_csv(a.related, header=None)[0].astype(str))
        drop(df[eid].astype(str).isin(rel), "related (user-supplied list)")

    drop(cov[[f"PC{i+1}" for i in range(a.n_pcs)]].isna().any(axis=1), "missing principal components")
    drop(cov.age.isna() | cov.sex.isna(), "missing age or sex")

    p = a.out_prefix
    cov.loc[keep].to_csv(f"{p}covariates.tsv", sep="\t", index=False, na_rep="NA")
    cov.loc[keep, ["FID", "IID"]].to_csv(f"{p}keep.tsv", sep="\t", index=False, header=False)

    txt = "\n".join([f"starting participants: {n0:,}", "exclusions:"] + log +
                    [f"\nanalysis sample: {int(keep.sum()):,}"])
    open(f"{p}exclusions.txt", "w").write(txt + "\n")
    print("\n" + txt)


if __name__ == "__main__":
    main()
