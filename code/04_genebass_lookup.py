#!/usr/bin/env python3
"""
00_genebass_lookup.py — look up candidate genes in Genebass (394,841 UK Biobank exomes).

Genebass has no documented public API, but the web app is backed by one:

    GET https://main.genebass.org/api/phenotypes
        -> 4,529 phenotypes with description, category, n_cases, analysis_id

    GET https://main.genebass.org/api/phewas/<ENSG>?burdenSet=pLoF
        -> {gene: ..., phewas: [ {phenocode, trait_type, coding, pheno_sex,
                                  Pvalue, Pvalue_Burden, Pvalue_SKAT,
                                  BETA_Burden, total_variants, ...}, ... ]}
        burdenSet is one of: pLoF | missense|LC | synonymous

RATE LIMITING MATTERS. Five concurrent requests got throttled hard, after which even
single requests hung for >140s. This script is deliberately serial with a delay, retries
with backoff, and caches every response so an interrupted run resumes without re-fetching.
Do not lower --delay below 1s. It is someone else's server.

Usage
    python3 00_genebass_lookup.py --genes ../genes/tier1.tsv --out genebass_renal.tsv
    python3 00_genebass_lookup.py --genes ../genes/ciliary.tsv --burden-set 'missense|LC'
"""
import argparse, gzip, io, json, os, sys, time, urllib.error, urllib.parse, urllib.request, zlib

API = "https://main.genebass.org/api"
UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 " \
     "(KHTML, like Gecko) Chrome/126.0 Safari/537.36"

def get(url, timeout=180, retries=4, delay=2.0):
    req = urllib.request.Request(url, headers={
        "User-Agent": UA, "Accept": "application/json",
        "Accept-Encoding": "gzip, deflate",
        "Referer": "https://app.genebass.org/", "Origin": "https://app.genebass.org",
    })
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                raw = r.read()
                enc = (r.headers.get("Content-Encoding") or "").lower()
                if enc == "gzip" or raw[:2] == b"\x1f\x8b":
                    raw = gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
                elif enc == "deflate":
                    try:
                        raw = zlib.decompress(raw)
                    except zlib.error:
                        raw = zlib.decompress(raw, -zlib.MAX_WBITS)
                return json.loads(raw.decode("utf-8"))
        except urllib.error.HTTPError as e:
            # 4xx other than 429 is permanent (gene absent, bad id) - do not retry
            if 400 <= e.code < 500 and e.code != 429:
                raise RuntimeError(f"HTTP {e.code} (permanent) for {url}") from e
            last = e
            wait = delay * (2 ** attempt)
            print(f"    retry {attempt+1}/{retries} in {wait:.0f}s (HTTP {e.code})",
                  file=sys.stderr, flush=True)
            time.sleep(wait)
            continue
        except Exception as e:
            last = e
            wait = delay * (2 ** attempt)
            print(f"    retry {attempt+1}/{retries} in {wait:.0f}s ({type(e).__name__})",
                  file=sys.stderr, flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed after {retries} attempts: {url}\n  {last}")


# Renal phenotype codes captured from /api/phenotypes (Genebass, 4,529 phenotypes).
# Hardcoded so the script does not depend on that endpoint, which is slow and times out.
# Matching is on (trait_type, phenocode) — phenocodes are unique per phenotype.
RENAL_PHENOS = {
    ("continuous", "30500"): ("Microalbumin in urine", 119013),
    ("continuous", "30510"): ("Creatinine (enzymatic) in urine", 383774),
    ("continuous", "30670"): ("Urea", 376551),
    ("continuous", "30700"): ("Creatinine", 376624),
    ("continuous", "30720"): ("Cystatin C", 376784),
    ("icd10", "C64"): ("C64 Malignant neoplasm of kidney", 985),
    ("icd_first_occurrence", "131290"): ("I12 hypertensive renal disease", 1513),
    ("icd_first_occurrence", "132006"): ("N04 nephrotic syndrome", 242),
    ("icd_first_occurrence", "132014"): ("N08 glomerular disorders in dis. class. elsewhere", 531),
    ("icd_first_occurrence", "132030"): ("N17 acute renal failure", 7312),
    ("icd_first_occurrence", "132032"): ("N18 chronic renal failure", 14086),
    ("icd_first_occurrence", "132034"): ("N19 unspecified renal failure", 2345),
    ("icd_first_occurrence", "132036"): ("N20 calculus of kidney and ureter", 5879),
    ("icd_first_occurrence", "132042"): ("N23 unspecified renal colic", 4304),
    ("icd_first_occurrence", "132044"): ("N25 impaired renal tubular function", 160),
    ("icd_first_occurrence", "132046"): ("N26 unspecified contracted kidney", 176),
    ("icd_first_occurrence", "132050"): ("N28 other disorders of kidney and ureter", 5275),
    ("icd_first_occurrence", "132530"): ("Q60 renal agenesis / reduction defects", 199),
    ("icd_first_occurrence", "132532"): ("Q61 cystic kidney disease", 735),
    ("icd_first_occurrence", "132536"): ("Q63 other congenital malformations of kidney", 452),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True,
                    help="TSV with gene_symbol and ensembl_id columns")
    ap.add_argument("--burden-set", default="pLoF",
                    choices=["pLoF", "missense|LC", "synonymous"])
    ap.add_argument("--delay", type=float, default=2.0, help="seconds between requests (min 1)")
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--timeout", type=float, default=180,
                    help="seconds per request; the API can be very slow")
    ap.add_argument("--mode", choices=["gene", "phenotype"], default="gene",
                    help="'gene': one request per gene (37-169 requests). "
                         "'phenotype': one request per renal trait (5 requests) that "
                         "returns ALL genes at once - far fewer requests, bigger payloads.")
    ap.add_argument("--all-phenotypes", action="store_true",
                    help="keep every phenotype, not just renal ones")
    ap.add_argument("--out", default="genebass_results.tsv")
    a = ap.parse_args()
    a.delay = max(a.delay, 1.0)
    os.makedirs(a.cache, exist_ok=True)

    genes = []
    with open(a.genes) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi, ei = hdr.index("gene_symbol"), hdr.index("ensembl_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(gi, ei) and f[ei].startswith("ENSG"):
                genes.append((f[gi], f[ei]))
    print(f"[genes] {len(genes)} from {a.genes}")

    if a.all_phenotypes:
        print("[pheno] fetching full phenotype dictionary (slow endpoint)")
        ph = get(f"{API}/phenotypes", timeout=a.timeout)
        ph = ph if isinstance(ph, list) else next(v for v in ph.values() if isinstance(v, list))
        keep = {(p.get("trait_type"), str(p.get("phenocode"))):
                (p.get("description"), p.get("n_cases")) for p in ph}
    else:
        keep = dict(RENAL_PHENOS)
    print(f"[pheno] {len(keep)} phenotypes retained")

    rows, failed = [], []
    bs = urllib.parse.quote(a.burden_set, safe="")
    want = {e: g for g, e in genes}
    want_sym = {g for g, _ in genes}

    if a.mode == "phenotype":
        # Continuous traits only: these carry the power, and the endpoint returns every
        # gene in one response. analysis_id = trait_type-phenocode-pheno_sex-coding-modifier
        CONT = [("30700", "Creatinine"), ("30720", "Cystatin C"), ("30670", "Urea"),
                ("30500", "Microalbumin in urine"),
                ("30510", "Creatinine (enzymatic) in urine")]
        for code, desc in CONT:
            aid = f"continuous-{code}-both_sexes--irnt"
            cf = os.path.join(a.cache, f"PHENO_{aid}_{a.burden_set.replace('|','_')}.json")
            if os.path.exists(cf):
                d = json.load(open(cf))
            else:
                print(f"[pheno-mode] {desc} ({aid})", flush=True)
                try:
                    d = get(f"{API}/phenotype/{aid}?burdenSet={bs}",
                            timeout=a.timeout, delay=a.delay)
                except Exception as e:
                    print(f"    FAILED: {e}", file=sys.stderr)
                    failed.append(desc)
                    continue
                json.dump(d, open(cf, "w"))
                time.sleep(a.delay)
            # defensive: locate the record array whatever the wrapper key is
            recs = d if isinstance(d, list) else next(
                (v for v in d.values() if isinstance(v, list)), [])
            n_hit = 0
            for r in recs:
                if not isinstance(r, dict):
                    continue
                sym = r.get("gene_symbol")
                eid = r.get("gene_id")
                if sym not in want_sym and eid not in want:
                    continue
                sym = sym or want.get(eid)
                n_hit += 1
                rows.append(dict(gene=sym, ensembl_id=eid or "", phenotype=desc,
                                 trait_type="continuous",
                                 n_cases=RENAL_PHENOS.get(("continuous", code), (desc, ""))[1],
                                 burden_set=a.burden_set,
                                 P_SKATO=r.get("Pvalue"), P_Burden=r.get("Pvalue_Burden"),
                                 P_SKAT=r.get("Pvalue_SKAT"), BETA=r.get("BETA_Burden"),
                                 n_variants=r.get("total_variants")))
            print(f"    {n_hit} of {len(genes)} target genes matched "
                  f"({len(recs):,} records in response)")
        genes = []   # skip the per-gene loop

    for i, (sym, ensg) in enumerate(genes, 1):
        cf = os.path.join(a.cache, f"{ensg}_{a.burden_set.replace('|','_')}.json")
        if os.path.exists(cf):
            d = json.load(open(cf))
        else:
            print(f"[{i:>4}/{len(genes)}] {sym} ({ensg})", flush=True)
            try:
                d = get(f"{API}/phewas/{ensg}?burdenSet={bs}", timeout=a.timeout, delay=a.delay)
            except Exception as e:
                print(f"    FAILED: {e}", file=sys.stderr)
                failed.append(sym)
                continue
            json.dump(d, open(cf, "w"))
            time.sleep(a.delay)
        for r in d.get("phewas", []):
            k = (r.get("trait_type"), str(r.get("phenocode")))
            if k in keep:
                desc, ncase = keep[k]
                rows.append(dict(gene=sym, ensembl_id=ensg, phenotype=desc,
                                 trait_type=r.get("trait_type"),
                                 n_cases=ncase, burden_set=a.burden_set,
                                 P_SKATO=r.get("Pvalue"), P_Burden=r.get("Pvalue_Burden"),
                                 P_SKAT=r.get("Pvalue_SKAT"), BETA=r.get("BETA_Burden"),
                                 n_variants=r.get("total_variants")))

    if not rows:
        sys.exit("No results. Check gene IDs and network access to main.genebass.org.")

    cols = ["gene", "ensembl_id", "phenotype", "trait_type", "n_cases", "burden_set",
            "P_SKATO", "P_Burden", "P_SKAT", "BETA", "n_variants"]
    rows.sort(key=lambda r: (r["P_SKATO"] is None, r["P_SKATO"]))
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    n_g = len({r["gene"] for r in rows})
    print(f"\n[write] {a.out}: {len(rows):,} rows, {n_g} genes")
    if failed:
        print(f"[warn] {len(failed)} genes failed: {', '.join(failed)}")
        print("       re-run to retry (cached genes are skipped)")

    # exome-wide significance in Genebass is ~2.5e-6 for gene-based tests
    sig = [r for r in rows if r["P_SKATO"] is not None and r["P_SKATO"] < 2.5e-6]
    print(f"\n[sig] {len(sig)} gene x renal-phenotype pairs at p < 2.5e-6:")
    for r in sig[:40]:
        print(f"  {r['gene']:<10} {str(r['phenotype'])[:44]:<44} "
              f"p={r['P_SKATO']:.2e}  beta={r['BETA']}")


if __name__ == "__main__":
    main()
