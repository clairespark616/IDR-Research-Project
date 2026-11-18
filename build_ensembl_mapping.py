"""
Build gene_symbol_to_ensembl.tsv directly from Ensembl's REST API
====================================================================
Replaces the UniProt-based approach, which can return transcript IDs
(ENST...) instead of gene IDs (ENSG...) depending on which column
was exported. This script queries Ensembl directly and only ever
returns gene-level IDs (ENSG...), which is what GTEx requires.

Usage:
  python3 build_ensembl_mapping.py

Inputs (same folder):
  sfari_cat1_genes.txt
  omim_cardiac_candidates.txt

Output:
  gene_symbol_to_ensembl.tsv   (overwrites any existing file)
"""

import requests
import time
import sys
import os

SFARI_FILE   = "sfari_cat1_genes.txt"
CARDIAC_FILE = "omim_cardiac_candidates.txt"
OUTPUT_FILE  = "gene_symbol_to_ensembl.tsv"

ENSEMBL_URL  = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}?content-type=application/json"
DELAY        = 0.2     # Ensembl allows ~15 requests/sec, this is safely under that
MAX_RETRIES  = 3
RETRY_WAIT   = 5

def load_genes(path):
    with open(path) as f:
        return [line.strip().upper() for line in f if line.strip()]

def fetch_ensembl_id(symbol):
    """
    Query Ensembl for the gene ID matching this symbol.
    Returns the ENSG... ID, or None if not found.
    """
    url = ENSEMBL_URL.format(symbol=symbol)
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=15)
            if r.status_code == 200:
                data = r.json()
                # Response is a list of cross-references; find the gene-type entry
                for entry in data:
                    if entry.get("id", "").startswith("ENSG"):
                        return entry["id"]
                return None
            elif r.status_code == 404:
                return None
            elif r.status_code == 429:
                time.sleep(RETRY_WAIT * attempt)
            else:
                time.sleep(RETRY_WAIT)
        except requests.exceptions.RequestException:
            time.sleep(RETRY_WAIT)
    return None

def main():
    print("\n" + "="*55)
    print("  Building gene_symbol_to_ensembl.tsv via Ensembl REST API")
    print("="*55)

    for f in [SFARI_FILE, CARDIAC_FILE]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    sfari   = load_genes(SFARI_FILE)
    cardiac = load_genes(CARDIAC_FILE)
    all_genes = sorted(set(sfari + cardiac))
    print(f"\nTotal unique genes to look up: {len(all_genes)}")

    # Resume support: skip genes already in the output file
    already_done = set()
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    already_done.add(parts[0].upper())
        print(f"Resuming — {len(already_done)} genes already done, skipping them.")

    to_fetch = [g for g in all_genes if g not in already_done]
    print(f"Genes to fetch: {len(to_fetch)}")

    if not to_fetch:
        print("All genes already fetched. Nothing to do.")
        return

    success = 0
    failed  = []

    with open(OUTPUT_FILE, "a") as out:
        for i, gene in enumerate(to_fetch, 1):
            ensg = fetch_ensembl_id(gene)
            if ensg:
                out.write(f"{gene}\t{ensg}\n")
                out.flush()
                print(f"[{i}/{len(to_fetch)}] {gene} -> {ensg}")
                success += 1
            else:
                print(f"[{i}/{len(to_fetch)}] {gene} -> NOT FOUND")
                failed.append(gene)
            time.sleep(DELAY)

    print("\n" + "="*55)
    print(f"Done. {success} genes mapped successfully.")
    if failed:
        print(f"{len(failed)} genes could not be mapped:")
        print(", ".join(failed))
        print("These genes will be excluded from the GTEx matching step.")
    print(f"Output written to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
