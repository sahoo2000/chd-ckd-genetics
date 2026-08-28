# Towards a Mechanistic Explanation of CHD-CKD

**Research Project 1 (BIOL61230)** — MSc Bioinformatics, The University of Manchester
Supervisor: Dr David Talavera

Do congenital heart disease (CHD) and chronic kidney disease (CKD) share genetic causes?

**Short answer: no.** Once you control for the fact that birth-defect syndromes damage many
organs at once, the apparent link disappears. What *does* share genetics with CHD is
congenital kidney *malformation*, and the mechanism is ciliary.

Everything here uses **public data only**.

> **Note on data versions.** This project has been re-analysed using the most recent releases of every resource available at the time of analysis (HPO 2026-06-23, gnomAD v4.1, MANE v1.5, FinnGen R10, MAGMA v1.10), rather than the versions current when the project was first conceived. Version numbers are stated throughout the thesis and in `docs/research_project_1.tex` Table 2 so that every result is exactly reproducible.


---

## Headline results

| Finding | Number |
|---|---|
| CHD to kidney malformation, adjusted for pleiotropy | OR 1.77 (1.52–2.07), p = 3.1e-13 |
| CHD to chronic kidney disease, adjusted for pleiotropy | **OR 0.95 (0.73–1.23), p = 0.68** |
| Genetic correlation, kidney function vs CHD | **rg = −0.03, p = 0.67** |
| Genetic correlation, kidney function vs CKD (positive control) | rg = −0.36, p = 6.1e-11 |
| Excess of rare-variant hits at p<0.001, ciliary genes | 9.7× |
| Direction of those hits | 59/60 toward worse kidney function, p = 5.3e-17 |

The full write-up is [`docs/research_project_1.pdf`](docs/research_project_1.pdf) — 33 pages, 16 figures, with declaration, acknowledgements, lists of figures and tables, a graphical abstract and five supplementary sections.

---

## Repository layout

```
chd-ckd-genetics/
├── README.md              <- you are here
├── code/                  <- run these in numerical order
│   ├── 01_download_data.sh
│   ├── 02_build_gene_list.py
│   ├── 03_hpo_analysis.py          <- the main analysis
│   ├── 04_genebass_lookup.py
│   ├── 05_genebass_analysis.py
│   ├── 06_ldsc.py
│   ├── 07_make_figures.py
│   ├── 08_run_magma.sh
│   ├── 09_make_schematic_figures.py
│   └── regenie_ukb/                <- unused, see note below
├── data/
│   ├── genes/             <- the candidate gene lists
│   └── results/           <- output tables
├── figures/               <- 16 PNG figures at 300 dpi
├── docs/
│   ├── research_project_1.tex
│   └── research_project_1.pdf
└── rawdata/               <- created by 01_download_data.sh (not in git)
```

---

## How to run it

You need Python 3 with `pandas`, `numpy`, `scipy`, `statsmodels` and `matplotlib`.

```bash
cd code
bash 01_download_data.sh          # ~5 GB, takes a while
python3 02_build_gene_list.py     # a few minutes
python3 03_hpo_analysis.py        # ~10 minutes
python3 04_genebass_lookup.py --genes ../data/genes/ciliary.tsv --out ciliary_renal.tsv
python3 05_genebass_analysis.py --results ciliary_renal.tsv
python3 06_ldsc.py --M 1173569 --trait "eGFR:egfr.hm3.tsv" ...
bash 08_run_magma.sh              # SLOW - about 4 hours
python3 07_make_figures.py
```

Then build the thesis:

```bash
cd docs
latexmk -pdf research_project_1.tex
```

**Timing warning.** `08_run_magma.sh` step 2 took about 4 hours on a laptop. Everything
else finishes in under an hour, except the downloads.

---

## What each script does

**`01_download_data.sh`** — fetches HPO, HGNC, MANE, gnomAD, CKDGen, FinnGen R10, the 1000
Genomes reference panel, MAGMA and the LD scores. Skips anything already downloaded.

**`02_build_gene_list.py`** — finds genes causing both a heart and a kidney malformation,
then cleans the list three ways: drops non-coding and non-autosomal genes, drops genes whose
only diseases kill in infancy (those people cannot be in an adult biobank), and separates out
immune genes that enter through autoimmune disease.

**`03_hpo_analysis.py`** — the main analysis, five tests:

1. Do heart and kidney genes overlap? (yes, OR 5.81)
2. Is the kidney special, or does every organ overlap the heart? (**it is not special** —
   kidney ranks 16th of 194 organ systems)
3. Why? Heart syndromes affect a median of 11 organ systems vs 5 for other diseases
4. Does the link survive adjusting for that? (kidney malformation yes, CKD no)
5. Which cardiac lesion carries the signal? (the severe ones, which adult cohorts lack)

**`04_genebass_lookup.py`** — queries [Genebass](https://app.genebass.org/) for rare-variant
burden results in 394,841 UK Biobank exomes. Genebass has no documented public API but the
web app is backed by one; the endpoints are recorded at the top of the script. **Please keep
the 2-second delay** — a concurrent burst got the API throttled during development.

**`05_genebass_analysis.py`** — tests whether signal remains after removing genes that are
already known, and whether the direction of effect is consistent.

**`06_ldsc.py`** — cross-trait LD score regression. The reference implementation is Python 2
and no longer installs, so this reimplements the published estimator. **It was validated
against a published heritability estimate and a positive control before use, but it is not
the reference binary — say so if you publish this.**

**`07_make_figures.py`** — draws the 11 data figures.

**`09_make_schematic_figures.py`** — draws the graphical abstract and the introduction and
methods diagrams.

**`08_run_magma.sh`** — gene-level and gene-set association from the CKDGen GWAS.

---

---

## Honest limitations

These are stated fully in the Discussion of the thesis. The short version:

- **HPO is Mendelian.** It describes rare syndromic disease and says nothing about
  common-variant kidney biology. `UMOD`, `APOL1` and `SHROOM3` are absent by construction.
- **The CHD genetic correlation is underpowered.** The heritability z-score for the FinnGen
  CHD endpoint is 2.12, below the conventional threshold of 4. The bound excluding
  |rg| > 0.23 is valid; a precise point estimate is not.
- **LD score regression used a reimplementation**, not the reference software.
- **The analysis is disease-level, not patient-level.** It measures whether genetic causes
  co-occur, not how often the conditions co-occur in living patients.
- **The subcomplex gradient rests on ~8 genes per group.** Treat it as hypothesis-generating.
- **CKDGen summary statistics are genomic-control corrected**, which deflates test
  statistics (observed heritability intercept 0.894). Slope-based rg is unaffected.

---

## Data sources

| Resource | Version | Access |
|---|---|---|
| Human Phenotype Ontology | 2026-06-23 | Open |
| HGNC complete set | current | Open |
| MANE Select | v1.5 | Open |
| gnomAD constraint | v4.1 | Open |
| Genebass | 394,841 UKB exomes | Open |
| CKDGen eGFR (Stanzick 2021) | European ancestry | Direct download |
| FinnGen | Release 10 | Direct download, no form |
| 1000 Genomes EUR panel | Phase 3 | Open |



---
