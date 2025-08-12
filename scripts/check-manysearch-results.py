#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from collections import OrderedDict

GENUS_ALIASES = {
    "Candida": ["Candida", "Candidozyma"],
}

SPECIES_ALIASES = {
    # Full "Genus species" -> aliases
    "Candida auris": ["Candida auris", "Candidozyma auris"],
}

def get_aliases(exp_genus: str, exp_species: str) -> tuple[list[str], list[str]]:
    """Return (genus_aliases, species_aliases). Species aliases include genus-swapped variants."""
    g_aliases = GENUS_ALIASES.get(exp_genus, [exp_genus])

    # Start from explicit species aliases if present, else the original full species
    s_aliases = set(SPECIES_ALIASES.get(exp_species, [exp_species]))

    # Add genus-swapped variants: keep epithet, swap genus through genus aliases
    if " " in exp_species:
        _, epithet = exp_species.split(" ", 1)
        for g in g_aliases:
            s_aliases.add(f"{g} {epithet}")

    return g_aliases, sorted(s_aliases)

def check_expected_in_text(text: str, exp_genus: str, exp_species: str) -> tuple[bool, bool]:
    """Return (genus_found, species_found) in `text`, allowing aliases."""
    if not isinstance(text, str):
        return False, False
    tl = text.lower()
    g_aliases, s_aliases = get_aliases(exp_genus, exp_species)
    genus_found = any(a.lower() in tl for a in g_aliases)
    species_found = any(a.lower() in tl for a in s_aliases)
    return genus_found, species_found


def prep_expected_taxonomy(expected_csv: str) -> dict:
    """
    Load expected taxonomy from CSV and prepare a mapping for quick lookup.
    """
    expected_df = pd.read_csv(expected_csv, header=None, names=["accession", "expected_genus", "expected_species"])
    # Replace underscores with spaces in species epithet
    expected_df["expected_species"] = expected_df["expected_species"].str.replace("_", " ", regex=False)
    expected_df["expected_species_full"] = expected_df["expected_genus"] + " " + expected_df["expected_species"]

    return OrderedDict(
        (row["accession"], {
            "genus": row["expected_genus"],
            "species": row["expected_species_full"]
        })
        for _, row in expected_df.iterrows()
    )

def load_bins_file(bins_file: str) -> dict:
    """
    Load bins file and return a mapping of accession to number of bins.
    """
    bins_per_acc = {}
    with open(bins_file) as f:
        for line in f:
            path = line.strip()
            if not path:
                continue
            acc = path.split("/")[2]  # e.g., SRR11734785 from multi/bins/SRR11734785/file.fa
            bins_per_acc[acc] = bins_per_acc.get(acc, 0) + 1
    return bins_per_acc



def main(args):

    # Read expected species list
    expected_map = prep_expected_taxonomy(args.expected_csv)
    
    # Load bins file list
    bins_per_acc = load_bins_file(args.bins)

    # Read raw metagenome manysearch results
    results_df = pd.read_csv(args.mgx_manysearch_csv)
   
    # Add columns for bin information
    results_df["total_n_bins"] = results_df["match_name"].map(lambda acc: bins_per_acc.get(acc, 0))
    results_df["no_bins"] = results_df["total_n_bins"] == 0

    # Add expected species and genus columns
    results_df["expected_species"] = None
    results_df["mgx_expected_species_found"] = False
    results_df["mgx_expected_genus_found"] = False

     # === Initialize bin columns ===

    # search bins x exact query genomes
    results_df["exact_bin_species_match"] = False
    results_df["exact_bin_genus_match"] = False
    results_df["exact_bin_best_ani"] = None
    results_df["exact_n_bins_with_match"] = 0

    # search bins x NCBI database
    results_df["ncbi_bin_species_match"] = False
    results_df["ncbi_bin_genus_match"] = False
    results_df["ncbi_bin_best_ani"] = None
    results_df["ncbi_n_bins_with_match"] = 0


    # Read bin manysearch results
    bin_manysearch_df = pd.read_csv(args.bin_manysearch_ncbi_csv)
    # Extract metagenome accession from bin query_name (SRRxxxx.y -> SRRxxxx)
    bin_manysearch_df["metagenome_accession"] = bin_manysearch_df["query_name"].astype(str).str.split(".", n=1).str[0]
    
    # read in bin multisearch
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

            # first, check exact genomes x metagenomes manysearch results
            # check if expected species and genus are found in the query name
            g_found, s_found = check_expected_in_text(row.get("query_name"), exp_genus, exp_species)
            results_df.at[idx, "mgx_expected_genus_found"] = g_found
            results_df.at[idx, "mgx_expected_species_found"] = s_found

            # --- Check bin-level matches ---
            if row["no_bins"]:
                # If no bins, skip bin checks
                continue

            ## Now check the multisearch (exact genomes x bins) results
            bin_rows = bins_df[bins_df["metagenome_accession"] == accession]

            # check if expected species and genus are found in bin matches
            g_aliases, s_aliases = get_aliases(exp_genus, exp_species)
            species_matches = bin_rows[bin_rows["match_name"].apply(
                lambda x: isinstance(x, str) and any(a.lower() in x.lower() for a in s_aliases)
            )]
            genus_matches = bin_rows[bin_rows["match_name"].apply(
                lambda x: isinstance(x, str) and any(a.lower() in x.lower() for a in g_aliases)
            )]

            if not species_matches.empty:
                results_df.at[idx, "exact_bin_species_match"] = True
                results_df.at[idx, "exact_n_bins_with_match"] = len(species_matches)
                results_df.at[idx, "exact_bin_best_ani"] = species_matches["average_containment_ani"].max()

            if not genus_matches.empty:
                results_df.at[idx, "exact_bin_genus_match"] = True
                # If species match is False but genus match exists, best ANI is from genus matches
                if not results_df.at[idx, "exact_bin_species_match"]:
                    results_df.at[idx, "exact_bin_best_ani"] = genus_matches["average_containment_ani"].max()
            
            # --- Check NCBI bin matches ---
            sub_ncbi = bin_manysearch_df[bin_manysearch_df["metagenome_accession"] == accession]
            if not sub_ncbi.empty:
                g_aliases, s_aliases = get_aliases(exp_genus, exp_species)

                s_rows_ncbi = sub_ncbi["match_name"].apply(
                    lambda x: isinstance(x, str) and any(a.lower() in x.lower() for a in s_aliases)
                )
                g_rows_ncbi = sub_ncbi["match_name"].apply(
                    lambda x: isinstance(x, str) and any(a.lower() in x.lower() for a in g_aliases)
                )

                species_rows_ncbi = sub_ncbi[s_rows_ncbi]
                genus_rows_ncbi = sub_ncbi[g_rows_ncbi]

                if not species_rows_ncbi.empty:
                    results_df.at[idx, "ncbi_bin_species_match"] = True
                    results_df.at[idx, "ncbi_n_bins_with_match"] = len(species_rows_ncbi)

                    # Prefer ANI if present; NCBI file has 'query_containment_ani'
                    ani_col = (
                        "query_containment_ani" if "query_containment_ani" in species_rows_ncbi
                        else ("containment" if "containment" in species_rows_ncbi else None)
                    )
                    if ani_col:
                        results_df.at[idx, "ncbi_bin_best_ani"] = species_rows_ncbi[ani_col].max()

                if not genus_rows_ncbi.empty:
                    results_df.at[idx, "ncbi_bin_genus_match"] = True
                    if not results_df.at[idx, "ncbi_bin_species_match"]:
                        ani_col = (
                            "query_containment_ani" if "query_containment_ani" in genus_rows_ncbi
                            else ("containment" if "containment" in genus_rows_ncbi else None)
                        )
                        if ani_col:
                            results_df.at[idx, "ncbi_bin_best_ani"] = genus_rows_ncbi[ani_col].max()


    # Sort results by expected list order
    acc_order = {acc: i for i, acc in enumerate(expected_map.keys())}
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
    wrong_species = results_df[(results_df["expected_species"].notna()) & (~results_df["mgx_expected_species_found"])]
    if not wrong_species.empty:
        print("\nAccessions with wrong species match:")
        for _, row in wrong_species.iterrows():
            print(f"  {row['match_name']} — expected {row['expected_species']}, got query {row['query_name']}")

    # Wrong genus
    wrong_genus = results_df[(results_df["expected_species"].notna()) & (~results_df["mgx_expected_genus_found"])]
    if not wrong_genus.empty:
        print("\nAccessions with wrong genus match:")
        for _, row in wrong_genus.iterrows():
            exp_genus = expected_map.get(row['match_name'], {}).get("genus", "N/A")
            print(f"  {row['match_name']} — expected genus {exp_genus}, got query {row['query_name']}")
    
    total_expected = len(expected_map)
    correct_species_count = results_df["mgx_expected_species_found"].sum()
    correct_genus_count = results_df["mgx_expected_genus_found"].sum()

    print("\n=== EXACT GENOMES X METAGENOMES ===")
    print(f"Total expected metagenome matches: {total_expected}")
    print(f"Correct species matches: {correct_species_count} / {total_expected}")
    print(f"Correct genus matches: {correct_genus_count} / {total_expected}")

    print("\n=== BINS X EXACT GENOMES ===")
    valid_bin_rows = results_df[~results_df["no_bins"]]
    bin_species_count = valid_bin_rows["exact_bin_species_match"].sum()
    bin_genus_count = valid_bin_rows["exact_bin_genus_match"].sum()
    print(f"Metagenomes with bin species match: {bin_species_count} / {len(valid_bin_rows)}")
    print(f"Metagenomes with bin genus match: {bin_genus_count} / {len(valid_bin_rows)}")

    print("\n=== BINS X NCBI GENOMES ===")
    ncbi_species_count = valid_bin_rows["ncbi_bin_species_match"].sum()
    ncbi_genus_count = valid_bin_rows["ncbi_bin_genus_match"].sum()
    print(f"Metagenomes with bin species match: {ncbi_species_count} / {len(valid_bin_rows)}")
    print(f"Metagenomes with bin genus match: {ncbi_genus_count} / {len(valid_bin_rows)}")

    missing_ncbi = valid_bin_rows[
        (~valid_bin_rows["ncbi_bin_species_match"]) & (~valid_bin_rows["ncbi_bin_genus_match"])
    ]
    if not missing_ncbi.empty:
        print("\nMetagenomes without any bin species/genus match (NCBI manysearch):")
        for _, r in missing_ncbi.iterrows():
            print(f"  {r['match_name']} — expected {expected_map[r['match_name']]['species']} ({r['total_n_bins']} bins)")


    # === Report missing bin matches ===
    missing_bin_matches = valid_bin_rows[
        (~valid_bin_rows["exact_bin_species_match"]) & (~valid_bin_rows["exact_bin_genus_match"])
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
    parser.add_argument("--mgx-manysearch-csv", help="CSV file with manysearch results.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.csv")
    parser.add_argument("--bin-manysearch-ncbi-csv", help="CSV file with bin search results.", default="output.manysearch/bins-x-ncbi-entire.manysearch.csv")
    parser.add_argument("--expected-csv", help="CSV file with expected metagenome accessions and species", default="multi/multi_mapping.csv")
    parser.add_argument("-o", "--output", help="Output CSV with added expected species/genus columns.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.checked.csv")
    parser.add_argument("--bins", help="path to all bins", default="multi_bins.txt")
    parser.add_argument("--multisearch-bins", help = "CSV file with multisearch bins", default="output.manysearch/bins-x-search-genomes.multisearch.default-thresh.csv")
    args = parser.parse_args()

    main(args)
