# 04_genebass_lookup.py
#
# Downloads rare variant results for our candidate genes from Genebass.
#
# Genebass has already done the hard work: it tested every gene against
# 4,529 different traits in 394,841 people who had their exome sequenced
# as part of UK Biobank. So we do not need to run the test ourselves, we
# just need to look up our genes.
#
# The Genebass website does not advertise a way to download data, but the
# website itself has to get the numbers from somewhere, and it uses these
# two web addresses:
#
#   https://main.genebass.org/api/phenotypes
#       gives the list of all 4,529 traits
#
#   https://main.genebass.org/api/phewas/<GENE>?burdenSet=pLoF
#       gives every result for one gene
#
# BE POLITE. This is a free academic server. Ask for one gene at a time
# and wait a couple of seconds in between. Sending lots of requests at
# once got the server to stop responding for several minutes while this
# script was being written.
#
# To run:   python3 04_genebass_lookup.py


import os
import time
import json
import gzip
import io
import urllib.request
import urllib.error


# ======================================================================
# SETTINGS
# ======================================================================

# which gene list to look up
GENE_FILE = "../data/genes/ciliary.tsv"

# where to save the answers
OUTPUT_FILE = "ciliary_renal.tsv"

# which kind of variants to ask about
#   "pLoF"        = variants that break the gene completely
#   "missense|LC" = variants that change one amino acid
BURDEN_SET = "pLoF"

# how long to wait between requests, in seconds. Do not go below 1.
SECONDS_BETWEEN_REQUESTS = 2

# how long to wait for one answer before giving up, in seconds
TIMEOUT = 300

# how many times to try again if a request fails
MAX_TRIES = 4

# downloaded answers are saved here so that if the script stops halfway
# through, running it again does not download the same genes twice
CACHE_FOLDER = "cache"

WEB_ADDRESS = "https://main.genebass.org/api"


# The kidney traits we care about. Genebass labels each trait with a type
# and a code, so we list them here rather than downloading the whole list
# of 4,529 traits every time. The numbers are how many people have each
# condition in UK Biobank.
KIDNEY_TRAITS = {
    ("continuous", "30500"): ("Microalbumin in urine", 119013),
    ("continuous", "30510"): ("Creatinine (enzymatic) in urine", 383774),
    ("continuous", "30670"): ("Urea", 376551),
    ("continuous", "30700"): ("Creatinine", 376624),
    ("continuous", "30720"): ("Cystatin C", 376784),
    ("icd10", "C64"): ("C64 Malignant neoplasm of kidney", 985),
    ("icd_first_occurrence", "131290"): ("I12 hypertensive renal disease", 1513),
    ("icd_first_occurrence", "132006"): ("N04 nephrotic syndrome", 242),
    ("icd_first_occurrence", "132014"): ("N08 glomerular disorders", 531),
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


# ======================================================================
# ASKING THE SERVER FOR ONE GENE
# ======================================================================

def download_one_gene(ensembl_id):
    """
    Ask Genebass for every result belonging to one gene.
    Returns the answer as a Python dictionary.
    """
    address = WEB_ADDRESS + "/phewas/" + ensembl_id + "?burdenSet=" + BURDEN_SET

    # Pretend to be a normal web browser, otherwise the server may refuse.
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/json",
        "Accept-Encoding": "gzip",
        "Referer": "https://app.genebass.org/",
    }
    request = urllib.request.Request(address, headers=headers)

    tries = 0
    while tries < MAX_TRIES:
        tries = tries + 1
        try:
            connection = urllib.request.urlopen(request, timeout=TIMEOUT)
            raw_bytes = connection.read()
            connection.close()

            # The server usually sends the answer compressed, so we have
            # to unzip it before we can read it. Compressed data always
            # starts with these two particular bytes.
            if raw_bytes[0:2] == b"\x1f\x8b":
                raw_bytes = gzip.GzipFile(fileobj=io.BytesIO(raw_bytes)).read()

            return json.loads(raw_bytes.decode("utf-8"))

        except urllib.error.HTTPError as error:
            # A 404 or similar means this gene simply is not there, so
            # there is no point trying again.
            if 400 <= error.code < 500 and error.code != 429:
                print("      gene not available (HTTP %d)" % error.code)
                return None
            wait = SECONDS_BETWEEN_REQUESTS * (2 ** tries)
            print("      server error %d, waiting %d seconds" % (error.code, wait))
            time.sleep(wait)

        except Exception:
            wait = SECONDS_BETWEEN_REQUESTS * (2 ** tries)
            print("      no answer, waiting %d seconds and trying again" % wait)
            time.sleep(wait)

    return None


# ======================================================================
# MAIN
# ======================================================================

# ---- read the gene list ----
gene_names = []
gene_ids = []

gene_file = open(GENE_FILE)
header_line = gene_file.readline().rstrip("\n").split("\t")
name_column = header_line.index("gene_symbol")
id_column = header_line.index("ensembl_id")

for line in gene_file:
    parts = line.rstrip("\n").split("\t")
    if len(parts) > id_column and parts[id_column].startswith("ENSG"):
        gene_names.append(parts[name_column])
        gene_ids.append(parts[id_column])
gene_file.close()

print("Looking up", len(gene_names), "genes from", GENE_FILE)
print("Variant type:", BURDEN_SET)
print("Kidney traits to keep:", len(KIDNEY_TRAITS))
print("")

if not os.path.exists(CACHE_FOLDER):
    os.makedirs(CACHE_FOLDER)

# ---- go through the genes one at a time ----
all_results = []
failed_genes = []

for position in range(len(gene_names)):
    gene_name = gene_names[position]
    gene_id = gene_ids[position]

    # Have we already downloaded this one?
    safe_name = BURDEN_SET.replace("|", "_")
    cache_file = os.path.join(CACHE_FOLDER, gene_id + "_" + safe_name + ".json")

    if os.path.exists(cache_file):
        answer = json.load(open(cache_file))
    else:
        print("  [%3d of %3d]  %s" % (position + 1, len(gene_names), gene_name))
        answer = download_one_gene(gene_id)

        if answer is None:
            failed_genes.append(gene_name)
            continue

        json.dump(answer, open(cache_file, "w"))
        time.sleep(SECONDS_BETWEEN_REQUESTS)

    # ---- keep only the kidney traits ----
    for one_result in answer.get("phewas", []):
        trait_type = one_result.get("trait_type")
        trait_code = str(one_result.get("phenocode"))

        if (trait_type, trait_code) in KIDNEY_TRAITS:
            trait_name, number_of_people = KIDNEY_TRAITS[(trait_type, trait_code)]

            all_results.append({
                "gene": gene_name,
                "ensembl_id": gene_id,
                "phenotype": trait_name,
                "trait_type": trait_type,
                "n_cases": number_of_people,
                "burden_set": BURDEN_SET,
                "P_SKATO": one_result.get("Pvalue"),
                "P_Burden": one_result.get("Pvalue_Burden"),
                "P_SKAT": one_result.get("Pvalue_SKAT"),
                "BETA": one_result.get("BETA_Burden"),
                "n_variants": one_result.get("total_variants"),
            })

# ---- write everything out ----
if len(all_results) == 0:
    print("")
    print("No results at all. Check the gene list and the internet connection.")
else:
    column_order = ["gene", "ensembl_id", "phenotype", "trait_type", "n_cases",
                    "burden_set", "P_SKATO", "P_Burden", "P_SKAT", "BETA",
                    "n_variants"]

    output = open(OUTPUT_FILE, "w")
    output.write("\t".join(column_order) + "\n")

    for result in all_results:
        row_pieces = []
        for column in column_order:
            value = result[column]
            if value is None:
                row_pieces.append("")
            else:
                row_pieces.append(str(value))
        output.write("\t".join(row_pieces) + "\n")
    output.close()

    genes_found = set()
    for result in all_results:
        genes_found.add(result["gene"])

    print("")
    print("Saved", len(all_results), "results for", len(genes_found),
          "genes to", OUTPUT_FILE)

if len(failed_genes) > 0:
    print("")
    print("These genes failed:", ", ".join(failed_genes))
    print("Run the script again to retry them. Genes already downloaded")
    print("are read from the cache folder, so nothing is fetched twice.")
