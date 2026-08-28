#!/usr/bin/env python3
"""
03_build_regenie_sets.py — emit regenie's --anno-file / --set-list / --mask-def.

Two set definitions are written from the same annotation:

  per-gene   one set per gene           -> N tests, Bonferroni alpha = 0.05/N
  pooled     ALL genes as a single set  -> 1 test, alpha = 0.05

The pooled set is the point. Collapsing 169 genes into one burden unit raises the
carrier frequency into a testable range AND drops the multiple-testing correction by
more than two orders of magnitude. Power to detect a 0.08 SD shift in eGFR, versus
0.21 SD per-gene.

Input: a VEP (or equivalent) annotation table with at least
    variant_id, gene_symbol, consequence
and optionally  revel, cadd_phred, loftee
Column names are configurable below.

regenie file formats (v3.x):
  anno file  : <variant_id> <set_name> <annotation>
  set list   : <set_name> <chr> <pos> <comma-separated variant ids>
  mask def   : <mask_name> <annotation>[,<annotation>...]
"""
import argparse, sys, os
from collections import defaultdict
import pandas as pd

# consequence terms treated as high-confidence loss of function
LOF_TERMS = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "start_lost", "stop_lost",
}
MISSENSE_TERMS = {"missense_variant"}


def classify(row, cols, revel_strict, revel_lenient):
    csq = str(row[cols["consequence"]]).split("&")
    is_lof = any(c in LOF_TERMS for c in csq)
    if is_lof:
        if cols.get("loftee") and str(row.get(cols["loftee"], "")).upper() == "LC":
            return None                      # low-confidence LoF: drop
        return "LoF"
    if any(c in MISSENSE_TERMS for c in csq):
        r = row.get(cols["revel"]) if cols.get("revel") else None
        try:
            r = float(r)
        except (TypeError, ValueError):
            return None
        if r >= revel_strict:
            return "missense_strict"
        if r >= revel_lenient:
            return "missense_lenient"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation", required=True, help="VEP-style variant annotation table")
    ap.add_argument("--genes", required=True, help="gene table (tsv with gene_symbol, chrom, gene_start)")
    ap.add_argument("--set-name", default="CANDIDATE_SET", help="name for the pooled set")
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--col-variant", default="variant_id")
    ap.add_argument("--col-gene", default="gene_symbol")
    ap.add_argument("--col-csq", default="consequence")
    ap.add_argument("--col-revel", default="revel")
    ap.add_argument("--col-loftee", default="loftee")
    ap.add_argument("--revel-strict", type=float, default=0.70)
    ap.add_argument("--revel-lenient", type=float, default=0.50)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    genes = pd.read_csv(a.genes, sep="\t")
    want = set(genes.gene_symbol)
    gpos = genes.set_index("gene_symbol")[["chrom", "gene_start"]].to_dict("index")
    print(f"[genes] {len(want)} target genes from {a.genes}")

    ann = pd.read_csv(a.annotation, sep=a.sep, low_memory=False)
    cols = {"variant": a.col_variant, "gene": a.col_gene, "consequence": a.col_csq}
    for k, v in [("revel", a.col_revel), ("loftee", a.col_loftee)]:
        if v in ann.columns:
            cols[k] = v
        else:
            print(f"[warn] column '{v}' absent — {k} filtering disabled")

    missing = [cols[k] for k in ("variant", "gene", "consequence") if cols[k] not in ann.columns]
    if missing:
        sys.exit(f"ERROR: required columns not in annotation file: {missing}")

    ann = ann[ann[cols["gene"]].isin(want)]
    print(f"[annot] {len(ann):,} variant records in target genes")

    anno_rows, per_gene, pooled = [], defaultdict(list), []
    counts = defaultdict(int)
    for _, r in ann.iterrows():
        cls = classify(r, cols, a.revel_strict, a.revel_lenient)
        if cls is None:
            continue
        vid, gene = str(r[cols["variant"]]), str(r[cols["gene"]])
        counts[cls] += 1
        anno_rows.append((vid, gene, cls))           # per-gene set membership
        anno_rows.append((vid, a.set_name, cls))     # pooled set membership
        per_gene[gene].append(vid)
        pooled.append(vid)

    if not pooled:
        sys.exit("ERROR: no qualifying variants. Check consequence terms and column names.")

    print("[class] " + "  ".join(f"{k}={v:,}" for k, v in sorted(counts.items())))
    print(f"[class] genes with >=1 qualifying variant: {len(per_gene)}/{len(want)}")
    dropped = sorted(want - set(per_gene))
    if dropped:
        print(f"[class] {len(dropped)} genes contribute nothing: {', '.join(dropped[:15])}"
              + (" ..." if len(dropped) > 15 else ""))

    p = lambda f: os.path.join(a.outdir, f)

    with open(p("annotations.txt"), "w") as fh:
        for vid, s, cls in anno_rows:
            fh.write(f"{vid} {s} {cls}\n")

    with open(p("setlist_pergene.txt"), "w") as fh:
        for g, vs in sorted(per_gene.items()):
            info = gpos.get(g)
            if not info:
                continue
            fh.write(f"{g} {info['chrom']} {int(info['gene_start'])} {','.join(vs)}\n")

    # pooled: regenie needs one representative chr/pos for the set; use the first gene's
    anchor = genes.sort_values(["chrom", "gene_start"]).iloc[0]
    with open(p("setlist_pooled.txt"), "w") as fh:
        fh.write(f"{a.set_name} {anchor.chrom} {int(anchor.gene_start)} {','.join(pooled)}\n")

    with open(p("masks.txt"), "w") as fh:
        fh.write("M1_LoF LoF\n")
        fh.write("M2_LoF_missense_strict LoF,missense_strict\n")
        fh.write("M3_LoF_missense_lenient LoF,missense_strict,missense_lenient\n")

    print(f"\n[write] {p('annotations.txt')}      ({len(anno_rows):,} rows)")
    print(f"[write] {p('setlist_pergene.txt')}   ({len(per_gene)} sets)")
    print(f"[write] {p('setlist_pooled.txt')}    (1 set, {len(pooled):,} variants)")
    print(f"[write] {p('masks.txt')}             (3 masks)")
    print(f"\nBonferroni alpha — per-gene: {0.05/max(len(per_gene),1):.3g}   pooled: 0.05")


if __name__ == "__main__":
    main()
