"""
IUPred3 API Fetcher
===================
Queries the IUPred3 REST API for every protein in a UniProt lengths TSV,
and writes the results in the exact format the main analysis script expects.

Run this TWICE — once for SFARI, once for cardiac:

  python3 fetch_iupred3.py sfari_cat1_lengths.tsv   sfari_iupred3_raw.txt
  python3 fetch_iupred3.py cardiac_lengths.tsv       cardiac_iupred3_raw.txt

Requirements:
  pip3 install requests

Inputs:
  lengths_tsv   — UniProt ID mapping TSV with at least a column containing
                  UniProt accessions (Entry column) and gene names
  output_file   — where to write the raw IUPred3 results

The script is resumable: if it gets interrupted partway through, re-run it
with the same arguments and it will skip proteins already in the output file.
"""

import requests
import time
import sys
import os

# ── Settings ────────────────────────────────────────────────────────────────

IUPRED3_URL   = "http://iupred3.elte.hu/iupred3/long/{accession}.json"
DELAY_SECONDS = 1.5   # wait between requests — be polite to the server
MAX_RETRIES   = 3     # retry failed requests this many times
RETRY_WAIT    = 10    # seconds to wait before retrying a failed request

# Browser-style headers — some academic servers reject plain script requests
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    ),
    "Accept": "application/json, text/plain, */*",
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def load_accessions(tsv_path):
    """
    Read the UniProt ID mapping TSV and return a list of
    (gene_symbol, uniprot_accession) tuples.

    The TSV must have an 'Entry' column (UniProt accession)
    and a gene name column ('Gene Names', 'From', etc.).
    """
    import csv

    pairs = []
    with open(tsv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        # Find the gene symbol column
        symbol_col = None
        for candidate in ["Gene Names", "Gene names", "From", "Gene Symbol"]:
            if candidate in reader.fieldnames:
                symbol_col = candidate
                break
        if symbol_col is None:
            print(f"ERROR: Could not find gene symbol column in {tsv_path}")
            print(f"Columns found: {reader.fieldnames}")
            sys.exit(1)

        if "Entry" not in reader.fieldnames:
            print(f"ERROR: Could not find 'Entry' (UniProt accession) column in {tsv_path}")
            print(f"Columns found: {reader.fieldnames}")
            sys.exit(1)

        for row in reader:
            accession   = row["Entry"].strip()
            gene_symbol = row[symbol_col].strip().upper()
            if accession and gene_symbol:
                # Gene Names can contain multiple symbols separated by spaces
                # Take the first one as the primary symbol
                gene_symbol = gene_symbol.split()[0]
                pairs.append((gene_symbol, accession))

    # Deduplicate by accession (keep first occurrence)
    seen = set()
    unique_pairs = []
    for pair in pairs:
        if pair[1] not in seen:
            seen.add(pair[1])
            unique_pairs.append(pair)

    print(f"Loaded {len(unique_pairs)} unique proteins from {tsv_path}")
    return unique_pairs


def load_already_done(output_path):
    """
    If the output file already exists (from a previous interrupted run),
    read it and return the set of accessions already processed.
    """
    done = set()
    if not os.path.exists(output_path):
        return done

    with open(output_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                # Header format: >sp|ACCESSION|...
                parts = line.strip().split("|")
                if len(parts) >= 2:
                    done.add(parts[1].strip())

    if done:
        print(f"Resuming — {len(done)} proteins already in {output_path}, skipping them.")
    return done


def fetch_iupred3(accession):
    """
    Query the IUPred3 REST API for one UniProt accession.
    Returns a list of (position, amino_acid, score) tuples,
    or None if the request failed after all retries.
    """
    url = IUPRED3_URL.format(accession=accession)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(url, headers=HEADERS, timeout=30)

            if response.status_code == 200:
                data = response.json()

                # IUPred3 JSON format:
                # {
                #   "sequence": "MASD...",
                #   "iupred2": [0.123, 0.456, ...],
                #   "anchor2": [...],
                #   "header": ">sp|Q9BYB0|SHANK3_HUMAN ..."
                # }
                sequence = data.get("sequence", "")
                scores   = data.get("iupred2", [])
                header   = data.get("header", f">sp|{accession}|")

                if not sequence or not scores:
                    print(f"  WARNING: Empty result for {accession} — skipping")
                    return None

                if len(sequence) != len(scores):
                    print(f"  WARNING: Sequence/score length mismatch for {accession} — skipping")
                    return None

                residues = [
                    (i + 1, aa, score)
                    for i, (aa, score) in enumerate(zip(sequence, scores))
                ]
                return header, residues

            elif response.status_code == 404:
                print(f"  WARNING: Accession {accession} not found in IUPred3 (404) — skipping")
                return None

            elif response.status_code == 429:
                # Rate limited — wait longer and retry
                print(f"  Rate limited. Waiting {RETRY_WAIT * attempt}s before retry...")
                time.sleep(RETRY_WAIT * attempt)

            else:
                print(f"  WARNING: HTTP {response.status_code} for {accession} "
                      f"(attempt {attempt}/{MAX_RETRIES})")
                time.sleep(RETRY_WAIT)

        except requests.exceptions.Timeout:
            print(f"  Timeout for {accession} (attempt {attempt}/{MAX_RETRIES})")
            time.sleep(RETRY_WAIT)

        except requests.exceptions.ConnectionError:
            print(f"  Connection error for {accession} (attempt {attempt}/{MAX_RETRIES})")
            time.sleep(RETRY_WAIT)

        except Exception as e:
            print(f"  Unexpected error for {accession}: {e} (attempt {attempt}/{MAX_RETRIES})")
            time.sleep(RETRY_WAIT)

    print(f"  FAILED after {MAX_RETRIES} attempts: {accession}")
    return None


def write_result(output_path, header, residues):
    """
    Append one protein's results to the output file in the standard
    three-column format the main analysis script expects.
    """
    with open(output_path, "a", encoding="utf-8") as f:
        f.write(header.strip() + "\n")
        for pos, aa, score in residues:
            f.write(f"{pos}\t{aa}\t{score:.4f}\n")
        f.write("\n")   # blank line between proteins


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 fetch_iupred3.py <lengths_tsv> <output_file>")
        print()
        print("Examples:")
        print("  python3 fetch_iupred3.py sfari_cat1_lengths.tsv   sfari_iupred3_raw.txt")
        print("  python3 fetch_iupred3.py cardiac_lengths.tsv       cardiac_iupred3_raw.txt")
        sys.exit(1)

    tsv_path    = sys.argv[1]
    output_path = sys.argv[2]

    # Load the list of proteins to process
    pairs = load_accessions(tsv_path)

    # Check which proteins are already done (for resuming)
    already_done = load_already_done(output_path)

    # Filter to only proteins not yet processed
    to_process = [(sym, acc) for sym, acc in pairs if acc not in already_done]
    print(f"Proteins to fetch: {len(to_process)}")

    if not to_process:
        print("All proteins already fetched. Nothing to do.")
        sys.exit(0)

    # Estimate runtime
    estimated_minutes = (len(to_process) * DELAY_SECONDS) / 60
    print(f"Estimated time: ~{estimated_minutes:.0f} minutes "
          f"(at {DELAY_SECONDS}s per protein)")
    print()

    # Fetch and write results
    success_count = 0
    fail_count    = 0

    for i, (gene_symbol, accession) in enumerate(to_process, 1):
        print(f"[{i}/{len(to_process)}] {gene_symbol} ({accession})...", end=" ")

        result = fetch_iupred3(accession)

        if result is not None:
            header, residues = result
            write_result(output_path, header, residues)
            disorder_fraction = sum(1 for _, _, s in residues if s > 0.5) / len(residues)
            print(f"done. Length={len(residues)}, disorder fraction={disorder_fraction:.3f}")
            success_count += 1
        else:
            print("skipped.")
            fail_count += 1

        # Be polite to the server
        if i < len(to_process):
            time.sleep(DELAY_SECONDS)

    # Summary
    print()
    print("=" * 50)
    print(f"Finished fetching IUPred3 data.")
    print(f"  Successful: {success_count}")
    print(f"  Failed/skipped: {fail_count}")
    print(f"  Output written to: {output_path}")

    if fail_count > 0:
        print()
        print(f"  Note: {fail_count} proteins could not be fetched.")
        print("  These will be excluded from the analysis automatically.")
        print("  Record them in your lab notebook as exclusions.")


if __name__ == "__main__":
    main()
