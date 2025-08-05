#!/usr/bin/env python3
import argparse
import pandas as pd

def main(args):
    # Load manifest
    df = pd.read_csv(args.manifest, comment="#")
    present_accessions = set(df['name'].unique())
    print(f"[INFO] Found {len(present_accessions)} unique accessions in manifest.")

    # Load accession list
    with open(args.accessions) as f:
        metagenomes = [line.strip() for line in f if line.strip()]
    print(f"[INFO] Loaded {len(metagenomes)} accessions from list.")

    # Find missing
    missing = [acc for acc in metagenomes if acc not in present_accessions]

    # Report
    if missing:
        print(f"[WARNING] {len(missing)} accessions are missing from manifest:")
        for acc in missing:
            print(f"  {acc}")
    else:
        print("[OK] All listed accessions are present in the manifest.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check for missing metagenome accessions in a SOURMASH manifest."
    )
    parser.add_argument(
        "--manifest",
        help="Path to SOURMASH manifest CSV file.",
        default="br-magtest.wort.mf.csv",
    )
    parser.add_argument(
        "--accessions",
        help="Path to file containing one accession per line.",
        default="multi/multinames.txt",
    )
    args = parser.parse_args()
    main(args)
