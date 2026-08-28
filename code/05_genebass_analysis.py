#!/usr/bin/env python3
"""
01_analyse_genebass.py — analyse the output of 00_genebass_lookup.py.

Answers four questions the raw table does not:
  1. Which genes clear exome-wide significance (2.5e-6)?
  2. Is there residual signal once those genes are removed, or is the rest noise?
  3. Does the residual signal point the right way (LoF should WORSEN kidney function)?
  4. Which ciliary subcomplex / gene family carries it?

The direction test is the important one. Selecting on p < 0.05 and then asking about the
sign of beta is a valid test: under the null, selected hits are still 50/50 in direction.
An excess of positive betas is therefore evidence of real signal that per-gene p-values,
individually underpowered, cannot show.
"""
import argparse, sys
import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu

EXOME_WIDE = 2.5e-6

# Known cystogenic / ADPKD genes. Their kidney signal is established; including them in a
# pooled test guarantees a positive result that rediscovers ADPKD rather than finding
# anything new.
CYSTOGENIC = {"PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11", "ALG5", "PKHD1"}

SUBCOMPLEX = {
    **{g: "BBSome" for g in ["BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBIP1", "TTC8"]},
    **{g: "Chaperonin (BBS6/10/12)" for g in ["MKKS", "BBS10", "BBS12"]},
    **{g: "IFT-A" for g in ["IFT122", "IFT140", "WDR19", "WDR35", "IFT43"]},
    **{g: "IFT-B" for g in ["IFT172", "IFT27", "IFT74", "IFT80"]},
    **{g: "Transition zone / MKS" for g in
       ["MKS1", "TMEM67", "TMEM216", "TMEM231", "TMEM237", "CC2D2A", "RPGRIP1L", "AHI1",
        "B9D1", "B9D2", "TMEM138", "TCTN1", "TCTN2"]},
    **{g: "NPHP module" for g in ["NPHP1", "NPHP3", "NPHP4", "INVS", "NEK8", "ANKS6", "SDCCAG8"]},
    **{g: "Dynein-2" for g in ["DYNC2H1", "DYNC2LI1", "DYNC2I1", "DYNC2I2"]},
    **{g: "ADPKD / cystogenic" for g in ["PKD1", "PKD2", "GANAB", "DNAJB11", "ALG9", "ALG5"]},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="output of 00_genebass_lookup.py")
    ap.add_argument("--exclude", nargs="*", default=sorted(CYSTOGENIC),
                    help="genes to remove for the residual-signal analysis")
    ap.add_argument("--out")
    a = ap.parse_args()

    d = pd.read_csv(a.results, sep="\t")
    d["P"] = pd.to_numeric(d.P_SKATO, errors="coerce")
    d["B"] = pd.to_numeric(d.BETA, errors="coerce")
    d = d.dropna(subset=["P"])
    if d.empty:
        sys.exit("No usable p-values in the results file.")
    bonf = 0.05 / len(d)

    print(f"tests {len(d):,} | genes {d.gene.nunique()} | phenotypes {d.phenotype.nunique()}")
    print(f"Bonferroni {bonf:.2g} | exome-wide {EXOME_WIDE:g}\n")

    # ---------------------------------------------------------------- 1. top hits
    sig = d[d.P < EXOME_WIDE].sort_values("P")
    print(f"=== exome-wide significant: {len(sig)} tests, {sig.gene.nunique()} genes ===")
    if len(sig):
        print(sig[["gene", "phenotype", "P", "BETA"]].head(25)
              .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    # ---------------------------------------------------------------- 2. residual
    excl = set(a.exclude)
    r = d[~d.gene.isin(excl)]
    print(f"\n=== residual signal, excluding {len(excl & set(d.gene))} known genes "
          f"({len(r):,} tests, {r.gene.nunique()} genes) ===")
    print(f"  {'threshold':<12}{'observed':>10}{'expected':>10}{'ratio':>8}")
    for thr in (0.05, 0.01, 1e-3, 1e-4):
        obs, exp = int((r.P < thr).sum()), thr * len(r)
        print(f"  p<{thr:<10g}{obs:>10}{exp:>10.1f}{obs/max(exp,1e-9):>8.1f}x")

    # ---------------------------------------------------------------- 3. direction
    print("\n=== direction (positive beta = rare LoF worsens kidney measure) ===")
    for lab, s in [("all tests", d), ("nominal p<0.05", d[d.P < 0.05]),
                   ("nominal, excl. known", r[r.P < 0.05])]:
        v = s.B.dropna()
        if not len(v):
            continue
        pos, n = int((v > 0).sum()), len(v)
        p = binomtest(pos, n, 0.5, alternative="greater").pvalue
        print(f"  {lab:<22}{pos:>4}/{n:<4} positive ({100*pos/n:>3.0f}%)  sign-test p={p:.3g}")

    # ---------------------------------------------------------------- 4. families
    d = d.copy()
    d["subcomplex"] = d.gene.map(SUBCOMPLEX).fillna("Unassigned")
    if d.subcomplex.nunique() > 1:
        print("\n=== enrichment by gene family ===")
        rows = []
        for sub, s in d.groupby("subcomplex"):
            nom, exp = int((s.P < 0.05).sum()), 0.05 * len(s)
            hits = s[s.P < 0.05]
            rows.append(dict(family=sub, genes=s.gene.nunique(), tests=len(s),
                             min_p=s.P.min(), nominal=nom, expected=exp,
                             enrichment=nom / exp if exp else np.nan,
                             pct_positive=100 * (hits.B > 0).mean() if len(hits) else np.nan))
        f = pd.DataFrame(rows).sort_values("enrichment", ascending=False)
        print(f.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    # ---------------------------------------------------------------- 5. power check
    if "n_variants" in d.columns:
        per = d.groupby("gene").agg(nv=("n_variants", "max"),
                                    nom=("P", lambda s: (s < 0.05).sum()))
        hi, lo = per[per.nom >= 2], per[per.nom < 2]
        if len(hi) > 2 and len(lo) > 2:
            p = mannwhitneyu(hi.nv, lo.nv, alternative="greater").pvalue
            print(f"\n=== is signal just power? ===")
            print(f"  >=2 nominal hits: {len(hi):>3} genes, median {hi.nv.median():.0f} variants")
            print(f"  < 2 nominal hits: {len(lo):>3} genes, median {lo.nv.median():.0f} variants")
            print(f"  Mann-Whitney p={p:.3g}" +
                  ("  -> signal is NOT explained by variant count" if p > 0.05 else
                   "  -> better-powered genes do show more signal; interpret with care"))

    # ---------------------------------------------------------------- 6. phenotypes
    print("\n=== phenotypes carrying residual signal (known genes excluded) ===")
    pp = r.groupby("phenotype").agg(genes=("gene", "nunique"),
                                    nominal=("P", lambda s: (s < 0.05).sum()),
                                    min_p=("P", "min"))
    pp["expected"] = 0.05 * pp.genes
    pp["enrichment"] = pp.nominal / pp.expected
    print(pp[pp.nominal > 0].sort_values("enrichment", ascending=False)
          .to_string(float_format=lambda x: f"{x:.3g}"))

    if a.out:
        d.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
        print(f"\n[write] {a.out}")


if __name__ == "__main__":
    main()
