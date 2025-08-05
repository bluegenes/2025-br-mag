#!/usr/bin/env python3
import argparse
import pandas as pd

def main(args):
    
    # Read expected species list
    expected_df = pd.read_csv(args.expected_csv, header=None, names=["accession", "expected_genus", "expected_species"])

    # Replace underscores with spaces in species epithet
    expected_df["expected_species"] = expected_df["expected_species"].str.replace("_", " ", regex=False)

    expected_df["expected_species_full"] = expected_df["expected_genus"] + " " + expected_df["expected_species"]

    # Create mapping for lookup
    expected_map = {
        row["accession"]: {
            "genus": row["expected_genus"],
            "species": row["expected_species_full"]
        }
        for _, row in expected_df.iterrows()
    }

    # Read results CSV
    results_df = pd.read_csv(args.manysearch_csv)

    # === Load bins file list and count bins per accession ===
    bins_per_acc = {}
    with open(args.bins) as f:
        for line in f:
            path = line.strip()
            if not path:
                continue
            acc = path.split("/")[2]  # e.g., SRR11734785 from multi/bins/SRR11734785/file.fa
            bins_per_acc[acc] = bins_per_acc.get(acc, 0) + 1

    # Add columns for bins
    results_df["total_n_bins"] = results_df["match_name"].map(lambda acc: bins_per_acc.get(acc, 0))
    results_df["no_bins"] = results_df["total_n_bins"] == 0

    # Add expected columns
    results_df["expected_species"] = None
    results_df["expected_species_found"] = False
    results_df["expected_genus_found"] = False

     # === Initialize bin columns ===
    results_df["bin_species_match"] = False
    results_df["bin_genus_match"] = False
    results_df["bin_best_ani"] = None
    results_df["n_bins_with_match"] = 0

    # reqad in bin multisearch
    bins_df = pd.read_csv(args.multisearch_bins)
    # Extract metagenome accession from bin query_name (SRRxxxx.y -> SRRxxxx)
    bins_df["metagenome_accession"] = bins_df["query_name"].str.split(".", n=1).str[0]

    # Check manysearch (metagenome) matches against expected species
    for idx, row in results_df.iterrows():
        accession = row["match_name"]  # match_name is the metagenome accession
        if accession in expected_map:
            exp_genus = expected_map[accession]["genus"]
            exp_species = expected_map[accession]["species"]

            results_df.at[idx, "expected_species"] = exp_species

            # Check genus match
            if pd.notna(row["query_name"]) and exp_genus.lower() in row["query_name"].lower():
                results_df.at[idx, "expected_genus_found"] = True

            # Check species match
            if pd.notna(row["query_name"]) and exp_species.lower() in row["query_name"].lower():
                results_df.at[idx, "expected_species_found"] = True


            # --- Check bin-level matches ---
            if row["no_bins"]:
                # If no bins, skip bin checks
                continue
            bin_rows = bins_df[bins_df["metagenome_accession"] == accession]
            species_matches = bin_rows[bin_rows["match_name"].str.contains(exp_species, case=False, na=False)]
            genus_matches = bin_rows[bin_rows["match_name"].str.contains(exp_genus, case=False, na=False)]

            if not species_matches.empty:
                results_df.at[idx, "bin_species_match"] = True
                results_df.at[idx, "n_bins_with_match"] = len(species_matches)
                results_df.at[idx, "bin_best_ani"] = species_matches["average_containment_ani"].max()

            if not genus_matches.empty:
                results_df.at[idx, "bin_genus_match"] = True
                # If species match is False but genus match exists, best ANI is from genus matches
                if not results_df.at[idx, "bin_species_match"]:
                    results_df.at[idx, "bin_best_ani"] = genus_matches["average_containment_ani"].max()

    # Sort results by expected list order
    acc_order = {acc: i for i, acc in enumerate(expected_df["accession"])}
    results_df["_sort_order"] = results_df["match_name"].map(lambda x: acc_order.get(x, len(acc_order)))
    results_df.sort_values("_sort_order", inplace=True)
    results_df.drop(columns=["_sort_order"], inplace=True)

    # Save updated CSV
    results_df.to_csv(args.output, index=False)
    print(f"Updated results written to '{args.output}'")

    # === REPORTING ===

    # Missing accessions
    found_accessions = set(results_df["match_name"])
    missing = [(acc, expected_map[acc]["species"]) for acc in expected_map if acc not in found_accessions]
    if missing:
        print("\nMissing accessions:")
        for acc, sp in missing:
            print(f"  {acc} — expected {sp}")
    else:
        print("\nAll metagenome accessions had genome matches.")

    # Wrong species
    wrong_species = results_df[(results_df["expected_species"].notna()) & (~results_df["expected_species_found"])]
    if not wrong_species.empty:
        print("\nAccessions with wrong species match:")
        for _, row in wrong_species.iterrows():
            print(f"  {row['match_name']} — expected {row['expected_species']}, got query {row['query_name']}")

    # Wrong genus
    wrong_genus = results_df[(results_df["expected_species"].notna()) & (~results_df["expected_genus_found"])]
    if not wrong_genus.empty:
        print("\nAccessions with wrong genus match:")
        for _, row in wrong_genus.iterrows():
            exp_genus = expected_map.get(row['match_name'], {}).get("genus", "N/A")
            print(f"  {row['match_name']} — expected genus {exp_genus}, got query {row['query_name']}")
    
    total_expected = len(expected_df)
    correct_species_count = results_df["expected_species_found"].sum()
    correct_genus_count = results_df["expected_genus_found"].sum()

    print("\n=== RAW METAGENOME MATCHES ===")
    print(f"Total expected metagenome matches: {total_expected}")
    print(f"Correct species matches: {correct_species_count} / {total_expected}")
    print(f"Correct genus matches: {correct_genus_count} / {total_expected}")

    print("\n=== BIN MATCHES ===")
    valid_bin_rows = results_df[~results_df["no_bins"]]
    bin_species_count = valid_bin_rows["bin_species_match"].sum()
    bin_genus_count = valid_bin_rows["bin_genus_match"].sum()
    print(f"Metagenomes with bin species match: {bin_species_count} / {len(valid_bin_rows)}")
    print(f"Metagenomes with bin genus match: {bin_genus_count} / {len(valid_bin_rows)}")

    # === Report missing bin matches ===
    missing_bin_matches = valid_bin_rows[
        (~valid_bin_rows["bin_species_match"]) & (~valid_bin_rows["bin_genus_match"])
    ]
    if not missing_bin_matches.empty:
        print("\nMetagenomes without any bin species/genus match:")
        for _, row in missing_bin_matches.iterrows():
            acc = row["match_name"]
            n_bins = row["total_n_bins"]
            print(f"  {acc} — expected {expected_map[acc]['species']} ({n_bins} bins)")
            bin_info = bins_df[bins_df["metagenome_accession"] == acc]
            if not bin_info.empty:
                print(bin_info.to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check metagenome search results against expected species list.")
    parser.add_argument("--manysearch-csv", help="CSV file with manysearch results.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.csv")
    parser.add_argument("--expected-csv", help="CSV file with expected metagenome accessions and species", default="multi/multi_mapping.csv")
    parser.add_argument("-o", "--output", help="Output CSV with added expected species/genus columns.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.checked.csv")
    parser.add_argument("--bins", help="path to all bins", default="multi_bins.txt")
    parser.add_argument("--multisearch-bins", help = "CSV file with multisearch bins", default="output.manysearch/bins-x-search-genomes.multisearch.default-thresh.csv")
    args = parser.parse_args()

    main(args)
