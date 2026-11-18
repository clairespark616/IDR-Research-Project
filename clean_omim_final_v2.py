"""
OMIM Cardiac Gene List Cleaner — Final Version
================================================
Hardcoded for genemap2.txt column structure confirmed via peek_omim.py:
  [6]  Gene/Locus And Other Related Symbols  <- gene symbol column
  [12] Phenotypes                             <- phenotype column

Usage:  python clean_omim_final.py

Inputs:
  genemap2.txt         (from omim.org/downloads)
  sfari_cat1_genes.txt (your SFARI Category 1 list)

Outputs:
  omim_cardiac_candidates.txt
  omim_cardiac_filtering_report.txt
"""

import pandas as pd
import re
import sys
import os

# ── File names ────────────────────────────────────────────────────────────────
OMIM_FILE   = "genemap2.txt"
SFARI_FILE  = "sfari_cat1_genes.txt"
OUTPUT_FILE = "omim_cardiac_candidates.txt"

# ── Hardcoded column names (confirmed from your file) ─────────────────────────
SYMBOL_COL    = "Gene/Locus And Other Related Symbols"
PHENOTYPE_COL = "Phenotypes"

# ── Keywords ──────────────────────────────────────────────────────────────────
CARDIAC_KEYWORDS = [
    "cardiomyopathy", "arrhythmia", "aortic", "channelopathy",
    "valve", "heart failure", "atrial fibrillation",
    "ventricular tachycardia", "ventricular fibrillation",
    "long qt", "brugada", "sick sinus", "dilated", "hypertrophic",
    "restrictive cardiomyopathy", "left ventricular",
    # Added to expand the candidate pool (see lab notebook note):
    "congenital heart", "cardiac conduction", "noncompaction",
    "myocarditis", "pericardial", "cardiac defect",
    "atrioventricular", "supravalvular", "coronary artery",
    "short qt", "catecholaminergic polymorphic",
    "cardiac arrest", "sudden cardiac death", "mitral",
    "tricuspid", "pulmonary valve", "septal defect",
    "conotruncal", "situs", "heterotaxy",
]

# Word-boundary matched — will NOT catch substrings like cardiomyopathy
NEURO_KEYWORDS = [
    "intellectual disability", "intellectual developmental",
    "autism spectrum", "autistic",
    "epilepsy", "epileptic", "encephalopathy",
    "neurodevelopmental", "developmental delay",
    "cognitive impairment", "ataxia",
    "peripheral neuropathy", "neurodegeneration", "neurodegenerative",
    "seizure disorder", "dementia", "schizophrenia", "psychosis",
    "cerebellar atrophy", "leukoencephalopathy", "spinocerebellar",
    "muscular dystrophy", "spinal muscular",
]

BAD_PREFIXES = [r"^\?", r"^\[", r"^\{"]

# ── Helpers ───────────────────────────────────────────────────────────────────
def contains_cardiac(text):
    t = str(text).lower()
    return any(kw in t for kw in CARDIAC_KEYWORDS)

def contains_neuro(text):
    t = str(text).lower()
    for kw in NEURO_KEYWORDS:
        if re.search(re.escape(kw) + r'\b', t):
            return True
    return False

def has_bad_prefix(text):
    for pheno in str(text).split(";"):
        pheno = pheno.strip()
        for pat in BAD_PREFIXES:
            if re.match(pat, pheno):
                return True
    return False

def load_omim(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        raw_lines = f.readlines()
    raw_lines[0] = raw_lines[0].lstrip("\ufeff")

    header_line = None
    header_idx  = None
    for i, line in enumerate(raw_lines):
        if line.startswith("#"):
            candidate = line.lstrip("#").strip()
            if "\t" in candidate:
                header_line = candidate
                header_idx  = i
        else:
            if header_line is not None:
                break

    if header_line is None:
        print("ERROR: Cannot find header in genemap2.txt")
        sys.exit(1)

    header = header_line.split("\t")
    n_cols = len(header)
    print(f"  Found {n_cols} columns. Gene symbol col index: "
          f"{header.index(SYMBOL_COL) if SYMBOL_COL in header else 'NOT FOUND'}")
    print(f"  Phenotype col index: "
          f"{header.index(PHENOTYPE_COL) if PHENOTYPE_COL in header else 'NOT FOUND'}")

    # Verify our hardcoded columns exist
    for col in [SYMBOL_COL, PHENOTYPE_COL]:
        if col not in header:
            print(f"\nERROR: Expected column '{col}' not found.")
            print(f"Columns in file: {header}")
            sys.exit(1)

    data_rows = []
    for line in raw_lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split("\t")
        if len(parts) < n_cols:
            parts += [""] * (n_cols - len(parts))
        else:
            parts = parts[:n_cols]
        data_rows.append(parts)

    df = pd.DataFrame(data_rows, columns=header)
    print(f"  Loaded {len(df)} rows")
    return df

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*55)
    print("  OMIM Cardiac Gene List Cleaner — Final")
    print("="*55)

    for f in [OMIM_FILE, SFARI_FILE]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found")
            sys.exit(1)

    # Load SFARI
    print(f"\n[Step 1] Loading SFARI gene list...")
    with open(SFARI_FILE) as f:
        sfari = {line.strip().upper() for line in f if line.strip()}
    print(f"  Loaded {len(sfari)} SFARI genes")

    # Load OMIM
    print(f"\n[Step 2] Loading genemap2.txt...")
    df = load_omim(OMIM_FILE)
    initial = len(df)

    # Show sample to confirm correct column
    sample_syms = df[SYMBOL_COL].dropna().head(5).tolist()
    print(f"  Sample gene symbols from '{SYMBOL_COL}': {sample_syms}")

    # Step 3: Filter
    print(f"\n[Step 3] Applying filters...")

    df = df[df[PHENOTYPE_COL].notna() & (df[PHENOTYPE_COL].str.strip() != "")]
    after_pheno = len(df)
    print(f"  Removed {initial - after_pheno:>5} rows: no phenotype")

    df = df[df[SYMBOL_COL].notna() & (df[SYMBOL_COL].str.strip() != "")]
    after_sym = len(df)
    print(f"  Removed {after_pheno - after_sym:>5} rows: no gene symbol")

    df = df[df[PHENOTYPE_COL].str.contains(r'\(3\)', regex=True, na=False)]
    after_key3 = len(df)
    print(f"  Removed {after_sym - after_key3:>5} rows: no Phenotype Mapping Key (3)")

    mask_bad = df[PHENOTYPE_COL].apply(has_bad_prefix)
    df = df[~mask_bad]
    after_prefix = len(df)
    print(f"  Removed {after_key3 - after_prefix:>5} rows: bad prefix (? [ {{)")

    mask_cardiac = df[PHENOTYPE_COL].apply(contains_cardiac)
    df = df[mask_cardiac]
    after_cardiac = len(df)
    print(f"  Removed {after_prefix - after_cardiac:>5} rows: no cardiac keyword")
    print(f"  Kept    {after_cardiac:>5} rows matching cardiac keywords")

    mask_neuro = df[PHENOTYPE_COL].apply(contains_neuro)
    df = df[~mask_neuro]
    after_neuro = len(df)
    print(f"  Removed {after_cardiac - after_neuro:>5} rows: neurological keyword")

    # Extract primary symbol — take first token before comma or space
    df = df.copy()
    df["primary_symbol"] = (
        df[SYMBOL_COL]
        .str.split(r"[,\s]+")
        .str[0]
        .str.strip()
        .str.upper()
    )
    df = df[df["primary_symbol"].notna() & (df["primary_symbol"] != "")]

    before_dedup = len(df)
    df = df.drop_duplicates(subset=["primary_symbol"])
    after_dedup = len(df)
    print(f"  Removed {before_dedup - after_dedup:>5} rows: duplicate symbol")

    mask_sfari = df["primary_symbol"].isin(sfari)
    n_overlap  = mask_sfari.sum()
    df = df[~mask_sfari]
    print(f"  Removed {n_overlap:>5} genes: overlap with SFARI")

    # Step 4: Results
    print(f"\n[Step 4] Results:")
    print(f"  Starting rows:            {initial}")
    print(f"  Final cardiac gene count: {len(df)}")

    if len(df) < 100:
        print(f"\n  WARNING: Only {len(df)} genes. Consider adding cardiac keywords.")

    symbols = sorted(df["primary_symbol"].tolist())
    print(f"\n  Genes ({len(symbols)} total):")
    for i in range(0, len(symbols), 10):
        print("    " + "  ".join(symbols[i:i+10]))

    # Step 5: Save
    print(f"\n[Step 5] Saving...")
    df["primary_symbol"].to_csv(OUTPUT_FILE, index=False, header=False)
    print(f"  Saved: {OUTPUT_FILE}")

    report = "omim_cardiac_filtering_report.txt"
    with open(report, "w") as f:
        f.write("OMIM Cardiac Gene Filtering Report\n")
        f.write("="*45 + "\n\n")
        f.write(f"Gene symbol column used: {SYMBOL_COL}\n")
        f.write(f"Phenotype column used:   {PHENOTYPE_COL}\n\n")
        f.write(f"Starting rows:                     {initial}\n")
        f.write(f"After no-phenotype removal:        {after_pheno}\n")
        f.write(f"After no-symbol removal:           {after_sym}\n")
        f.write(f"After Key(3) filter:               {after_key3}\n")
        f.write(f"After bad-prefix removal:          {after_prefix}\n")
        f.write(f"After cardiac keyword filter:      {after_cardiac}\n")
        f.write(f"After neurological keyword filter: {after_neuro}\n")
        f.write(f"After deduplication:               {after_dedup}\n")
        f.write(f"After SFARI overlap removal:       {len(df)}\n\n")
        f.write(f"Final cardiac gene count:          {len(df)}\n\n")
        f.write("Cardiac keywords:\n")
        for kw in CARDIAC_KEYWORDS:
            f.write(f"  - {kw}\n")
        f.write("\nNeurological exclusion keywords:\n")
        for kw in NEURO_KEYWORDS:
            f.write(f"  - {kw}\n")
        f.write("\nFinal gene list:\n")
        for s in symbols:
            f.write(f"  {s}\n")

    print(f"  Saved: {report}\n")

if __name__ == "__main__":
    main()
