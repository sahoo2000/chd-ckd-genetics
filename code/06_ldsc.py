#!/usr/bin/env python3
"""
ldsc_lite.py — LD Score regression (heritability and cross-trait genetic correlation).

A self-contained implementation of Bulik-Sullivan et al. (Nat Genet 2015a, 2015b). The
reference implementation is Python 2 only, which no longer installs cleanly; this follows
the same estimators. It is NOT the reference binary — validate the h2 it produces against
a published estimate for at least one trait before trusting the rg.

Model
  h2:  E[chi2_j]    = 1 + (N h2 / M) * l_j
  rg:  E[z1j z2j]   = (sqrt(N1 N2) rho_g / M) * l_j + intercept

The intercept in the rg regression absorbs sample overlap, so two GWAS from the same
cohort can be compared without biasing rho_g.

Standard errors come from a 200-block genome-ordered jackknife. Weights follow the
two-step scheme: an LD weight 1/max(l,1) for the correlation induced between nearby SNPs,
and a heteroscedasticity weight from a first-pass h2 estimate.

Genetic correlation is invariant to the observed/liability scale transformation, so rg is
reported without any prevalence assumption. Observed-scale h2 for a binary trait is NOT
comparable across traits with different prevalence — treat it as an input to rg, not a
result.
"""
import argparse, sys
import numpy as np
import pandas as pd

AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
N_BLOCKS = 200


def munge(path, name, n=None, ncase=None, nctrl=None, maf_min=0.01):
    """Read a filtered GWAS: SNP A1 A2 BETA SE [N] [FREQ]. A1 is the effect allele."""
    cols = ["SNP", "A1", "A2", "BETA", "SE", "C6", "C7"]
    d = pd.read_csv(path, sep="\t", header=None, names=cols, engine="c",
                    dtype={"SNP": str, "A1": str, "A2": str})
    n0 = len(d)
    if n is not None:
        d["N"] = float(n)
        d["FRQ"] = pd.to_numeric(d.C6, errors="coerce")
    else:
        d["N"] = pd.to_numeric(d.C6, errors="coerce")
        d["FRQ"] = pd.to_numeric(d.C7, errors="coerce")
    if ncase is not None and nctrl is not None:
        d["N"] = 4.0 / (1.0 / ncase + 1.0 / nctrl)      # effective N for binary traits

    d["BETA"] = pd.to_numeric(d.BETA, errors="coerce")
    d["SE"] = pd.to_numeric(d.SE, errors="coerce")
    d = d.dropna(subset=["BETA", "SE", "N"])
    d = d[d.SE > 0]
    d["Z"] = d.BETA / d.SE

    d = d[~d.SNP.duplicated(keep=False)]                 # drop any duplicated rsID
    amb = [(a, b) in AMBIGUOUS for a, b in zip(d.A1, d.A2)]
    d = d[~np.array(amb)]                                # strand-ambiguous
    d = d[d.A1.isin(list("ACGT")) & d.A2.isin(list("ACGT"))]
    if d.FRQ.notna().any():
        d = d[(d.FRQ > maf_min) & (d.FRQ < 1 - maf_min)]

    chi2_max = max(80.0, 0.001 * d.N.max())
    n_out = int((d.Z ** 2 > chi2_max).sum())
    d = d[d.Z ** 2 <= chi2_max]
    print(f"  [{name}] {n0:,} -> {len(d):,} SNPs "
          f"(dropped {n_out} with chi2 > {chi2_max:.0f}); mean N = {d.N.mean():,.0f}")
    return d[["SNP", "A1", "A2", "Z", "N"]]


def _wls(x, y, w):
    X = np.column_stack([np.ones_like(x), x])
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    return beta  # [intercept, slope]


def _jackknife(x, y, w, fn):
    """200-block genome-ordered jackknife; returns (estimate, se)."""
    n = len(x)
    full = fn(_wls(x, y, w))
    idx = np.array_split(np.arange(n), min(N_BLOCKS, n))
    vals = []
    for b in idx:
        keep = np.ones(n, bool)
        keep[b] = False
        vals.append(fn(_wls(x[keep], y[keep], w[keep])))
    vals = np.array(vals, float)
    k = len(vals)
    pseudo = k * full - (k - 1) * vals
    return full, pseudo.std(ddof=1) / np.sqrt(k)


def estimate_h2(d, ld, M):
    """Observed-scale SNP heritability."""
    m = d.merge(ld, on="SNP")
    l, chi2, N = m.L2.values, m.Z.values ** 2, m.N.values
    lw = np.maximum(l, 1.0)
    x = l * N / M

    h2_init = max(min(_wls(x, chi2, 1.0 / lw)[1], 1.0), 0.0)     # step 1
    het = (1.0 + h2_init * x) ** 2                                # step 2 weights
    w = 1.0 / (lw * het)

    est, se = _jackknife(x, chi2, w, lambda b: b[1])
    icpt, icpt_se = _jackknife(x, chi2, w, lambda b: b[0])
    return dict(h2=est, se=se, intercept=icpt, intercept_se=icpt_se,
                mean_chi2=float(chi2.mean()), n_snp=len(m), mean_N=float(N.mean()))


def estimate_rg(d1, d2, ld, M, h2_1, h2_2):
    """Cross-trait genetic covariance and correlation, allele-aligned on trait 1."""
    m = d1.merge(d2, on="SNP", suffixes=("_1", "_2")).merge(ld, on="SNP")
    same = (m.A1_1 == m.A1_2) & (m.A2_1 == m.A2_2)
    flip = (m.A1_1 == m.A2_2) & (m.A2_1 == m.A1_2)
    m = m[same | flip].copy()
    m.loc[flip[same | flip].values, "Z_2"] *= -1     # align to trait 1's effect allele

    l = m.L2.values
    z1, z2 = m.Z_1.values, m.Z_2.values
    N1, N2 = m.N_1.values, m.N_2.values
    lw = np.maximum(l, 1.0)
    Nbar = np.sqrt(N1 * N2)
    x = l * Nbar / M

    rho_init = _wls(x, z1 * z2, 1.0 / lw)[1]
    v1 = h2_1 * l * N1 / M + 1.0
    v2 = h2_2 * l * N2 / M + 1.0
    w = 1.0 / (lw * (v1 * v2 + (rho_init * x) ** 2))

    rho, rho_se = _jackknife(x, z1 * z2, w, lambda b: b[1])
    icpt, icpt_se = _jackknife(x, z1 * z2, w, lambda b: b[0])
    denom = np.sqrt(max(h2_1, 1e-12) * max(h2_2, 1e-12))
    return dict(rho_g=rho, rho_g_se=rho_se, rg=rho / denom, rg_se=rho_se / denom,
                intercept=icpt, intercept_se=icpt_se, n_snp=len(m))


def load_spec(spec):
    """name:path[:ncase,nctrl] — omit counts for a continuous trait with an N column."""
    parts = spec.split(":")
    name, path = parts[0], parts[1]
    ncase = nctrl = None
    if len(parts) > 2 and parts[2]:
        ncase, nctrl = (float(v) for v in parts[2].split(","))
    return name, path, ncase, nctrl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", action="append", required=True,
                    help="name:path[:ncase,nctrl]; repeatable")
    ap.add_argument("--ldscores", default="ldscores.tsv.gz")
    ap.add_argument("--M", type=int, required=True, help="M_5_50 across autosomes")
    ap.add_argument("--out", default="ldsc_results.tsv")
    a = ap.parse_args()

    ld = pd.read_csv(a.ldscores, sep="\t")[["SNP", "L2"]]
    print(f"[ld] {len(ld):,} LD-score SNPs; M_5_50 = {a.M:,}\n[munge]")

    traits, h2 = {}, {}
    for spec in a.trait:
        name, path, nc, nk = load_spec(spec)
        traits[name] = munge(path, name, ncase=nc, nctrl=nk)

    print("\n=== SNP heritability (observed scale) ===")
    rows = []
    for name, d in traits.items():
        r = estimate_h2(d, ld, a.M)
        h2[name] = r["h2"]
        z = r["h2"] / r["se"] if r["se"] else np.nan
        print(f"  {name:<24} h2 = {r['h2']:.4f} ({r['se']:.4f})  z = {z:5.2f}   "
              f"intercept = {r['intercept']:.3f} ({r['intercept_se']:.3f})  "
              f"mean chi2 = {r['mean_chi2']:.3f}  SNPs = {r['n_snp']:,}")
        rows.append(dict(analysis="h2", trait1=name, trait2="", estimate=r["h2"],
                         se=r["se"], z=z, intercept=r["intercept"], n_snp=r["n_snp"]))
        if z < 4:
            print(f"      NOTE: h2 z-score below 4 — LDSC genetic correlation involving "
                  f"{name} is unreliable.")

    print("\n=== genetic correlation ===")
    names = list(traits)
    from scipy.stats import norm
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            r = estimate_rg(traits[n1], traits[n2], ld, a.M, h2[n1], h2[n2])
            z = r["rg"] / r["rg_se"] if r["rg_se"] else np.nan
            p = 2 * norm.sf(abs(z)) if np.isfinite(z) else np.nan
            print(f"  {n1} ~ {n2}")
            print(f"      rg = {r['rg']:.4f} ({r['rg_se']:.4f})  z = {z:.2f}  p = {p:.3g}"
                  f"   [{r['rg']-1.96*r['rg_se']:.3f}, {r['rg']+1.96*r['rg_se']:.3f}]")
            print(f"      rho_g = {r['rho_g']:.5f} ({r['rho_g_se']:.5f})   "
                  f"intercept = {r['intercept']:.4f} ({r['intercept_se']:.4f})  "
                  f"SNPs = {r['n_snp']:,}")
            rows.append(dict(analysis="rg", trait1=n1, trait2=n2, estimate=r["rg"],
                             se=r["rg_se"], z=z, intercept=r["intercept"],
                             n_snp=r["n_snp"]))

    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(f"\n[write] {a.out}")


if __name__ == "__main__":
    main()
