"""
Build gene_symbol_to_ensembl.tsv via mygene.info (batch lookup)
==================================================================
Replaces the one-gene-at-a-time Ensembl REST approach with a single
batch query to mygene.info, which searches across symbols, aliases,
and previous/withdrawn names — giving a much higher match rate.

Usage:  python3 build_ensembl_mygene.py

Inputs (same folder):
  sfari_cat1_genes.txt
  omim_cardiac_candidates.txt

Output:
  gene_symbol_to_ensembl.tsv
  unmapped_genes.txt   (genes mygene.info could not resolve, for manual review)
"""

import requests
import sys
import os
import json

SFARI_FILE   = "sfari_cat1_genes.txt"
CARDIAC_FILE = "omim_cardiac_candidates.txt"
OUTPUT_FILE  = "gene_symbol_to_ensembl.tsv"
UNMAPPED_FILE = "unmapped_genes.txt"

MYGENE_URL = "https://mygene.info/v3/query"
BATCH_SIZE = 100   # mygene.info handles large batches well; chunk for safety

def load_genes(path):
    with open(path) as f:
        return [line.strip().upper() for line in f if line.strip()]

def query_batch(symbols):
    """
    Query mygene.info for a batch of gene symbols.
    Searches symbol, alias, and previous/withdrawn names (scopes parameter).
    Returns a dict: symbol -> ENSG id (or None if not found)
    """
    payload = {
        "q": ",".join(symbols),
        "scopes": "symbol,alias,name,retired",
        "species": "human",
        "fields": "ensembl.gene,symbol,query",
    }
    try:
        r = requests.post(MYGENE_URL, data=payload, timeout=30)
        r.raise_for_status()
        results = r.json()
    except requests.exceptions.RequestException as e:
        print(f"  ERROR querying batch: {e}")
        return {s: None for s in symbols}

    mapping = {}
    for entry in results:
        query = entry.get("query", "").upper()
        if entry.get("notfound"):
            mapping[query] = None
            continue

        ensembl = entry.get("ensembl")
        ensg = None
        if isinstance(ensembl, dict):
            gene = ensembl.get("gene")
            if isinstance(gene, list):
                ensg = gene[0]
            else:
                ensg = gene
        elif isinstance(ensembl, list):
            # Multiple entries — take first one with a gene ID
            for e in ensembl:
                gene = e.get("gene")
                if gene:
                    ensg = gene[0] if isinstance(gene, list) else gene
                    break

        # Only keep first/best match per query symbol
        if query not in mapping or (mapping.get(query) is None and ensg):
            mapping[query] = ensg

    # Make sure every requested symbol has an entry (even if None)
    for s in symbols:
        if s not in mapping:
            mapping[s] = None

    return mapping

def main():
    print("\n" + "="*55)
    print("  Building gene_symbol_to_ensembl.tsv via mygene.info")
    print("="*55)

    for f in [SFARI_FILE, CARDIAC_FILE]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    sfari   = load_genes(SFARI_FILE)
    cardiac = load_genes(CARDIAC_FILE)
    all_genes = sorted(set(sfari + cardiac))
    print(f"\nTotal unique genes to look up: {len(all_genes)}")

    all_mappings = {}
    unmapped = []

    # Process in batches
    for i in range(0, len(all_genes), BATCH_SIZE):
        batch = all_genes[i:i + BATCH_SIZE]
        print(f"\nQuerying batch {i//BATCH_SIZE + 1} "
              f"({i+1}-{min(i+BATCH_SIZE, len(all_genes))} of {len(all_genes)})...")
        batch_result = query_batch(batch)
        all_mappings.update(batch_result)

    # Separate found vs unmapped
    found = {k: v for k, v in all_mappings.items() if v}
    unmapped = sorted([k for k, v in all_mappings.items() if not v])

    print(f"\n{'='*55}")
    print(f"Mapped:   {len(found)} / {len(all_genes)}")
    print(f"Unmapped: {len(unmapped)} / {len(all_genes)}")

    # Write output
    with open(OUTPUT_FILE, "w") as out:
        for symbol, ensg in sorted(found.items()):
            out.write(f"{symbol}\t{ensg}\n")
    print(f"\nSaved mappings to: {OUTPUT_FILE}")

    if unmapped:
        with open(UNMAPPED_FILE, "w") as out:
            for symbol in unmapped:
                out.write(f"{symbol}\n")
        print(f"Saved unmapped gene list to: {UNMAPPED_FILE}")
        print(f"\nUnmapped genes (first 30 shown):")
        print(", ".join(unmapped[:30]) + (" ..." if len(unmapped) > 30 else ""))
        print(f"\nThese genes will be excluded from GTEx matching.")
        print(f"Record this count in your lab notebook for the Methods section.")
        print(f"\nFor any critical genes you want to manually recover, search")
        print(f"the symbol on genenames.org (HGNC) to find its current Ensembl ID,")
        print(f"then add a line to {OUTPUT_FILE} in the format:  SYMBOL<tab>ENSG...")

if __name__ == "__main__":
    main()
