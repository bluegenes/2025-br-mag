#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from collections import OrderedDict

GENUS_ALIASES = {
    "Candida": ["Candida", "Candidozyma"],
    "Acanthamoeba": ["Mimivirus"], # NOT Acanthamoeba
}

SPECIES_ALIASES = {
    # Full "Genus species" -> aliases
    "Candida auris": ["Candida auris", "Candidozyma auris"],
    "Acanthamoeba polyphaga mimivirus": ["Acanthamoeba polyphaga mimivirus", "Mimivirus bradfordmassiliense"]
}

NCBI_EXCLUDE_GENUS = {g.strip().lower() for g in ["Acanthamoeba"]} 

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
    tl = text.lower().strip()
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
            "genus": row["expected_genus"].strip(),
            "species": row["expected_species_full"].strip()
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


# prefer ANI if present; fall back to containment
def get_top_match_info(df):
    if "query_containment_ani" in df:
        max_row = df.loc[df["query_containment_ani"].idxmax()]
        return max_row["query_containment_ani"], max_row["match_name"]
    if "average_containment_ani" in df:
        max_row = df.loc[df["average_containment_ani"].idxmax()]
        return max_row["average_containment_ani"], max_row["match_name"]
    if "containment" in df:
        max_row = df.loc[df["containment"].idxmax()]
        return max_row["containment"], max_row["match_name"]
    return None, None


def check_bat_for_accession(acc: str, exp_genus: str, exp_species: str, bat_dir="output.BAT"):
    """
    Return (bat_bin_match_level, bat_bin_support, bat_orf_match_level) for a given accession.
    """
    import os, re

    bin_match_level = "unmatched"
    bin_support = None
    orf_match_level = "unmatched"

    acc_dir = os.path.join(bat_dir, acc)
    bin_file = os.path.join(acc_dir, "out.BAT.bin2classification.taxnames.txt")
    orf_file = os.path.join(acc_dir, "out.BAT.ORF2LCA.taxnames.txt")

    # Expand aliases
    genus_aliases, species_aliases = get_aliases(exp_genus, exp_species)

    # --- Bin-level classifications ---
    acc = acc.strip()
    print(acc)
    if acc == "SRR3458563":
        import pdb;pdb.set_trace()

    if os.path.exists(bin_file):
        try:
            with open(bin_file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 5:
                        continue
                    classification_text = " ".join(parts[5:])
                    if acc == "SRR3458563":
                        import pdb;pdb.set_trace()

                    # --- check species aliases first ---
                    species_scores = []
                    for alias in species_aliases:
                        m = re.search(rf"{re.escape(alias)}:\s*([\d.]+)", classification_text, re.I)
                        if m:
                            species_scores.append(float(m.group(1)))
                    if species_scores:
                        # species always wins; take highest score
                        best_score = max(species_scores)
                        if bin_match_level != "species" or (bin_support is None or best_score > bin_support):
                            bin_match_level, bin_support = "species", best_score
                        continue  # no need to check genus for this line

                    # --- otherwise check genus aliases ---
                    genus_scores = []
                    for alias in genus_aliases:
                        m = re.search(rf"{re.escape(alias)}:\s*([\d.]+)", classification_text, re.I)
                        if m:
                            genus_scores.append(float(m.group(1)))
                    if genus_scores:
                        best_score = max(genus_scores)
                        if bin_match_level != "species":  # don't override species
                            if bin_match_level != "genus" or (bin_support is None or best_score > bin_support):
                                bin_match_level, bin_support = "genus", best_score
        except Exception as e:
            print(f"Warning: could not parse {bin_file}: {e}")
    else:
        bin_match_level = "NA"

    # --- ORF-level classifications ---
    if os.path.exists(orf_file):
        try:
            with open(orf_file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    # species aliases first
                    if any(alias.lower() in line.lower() for alias in species_aliases):
                        orf_match_level = "species"
                        break
                    # genus aliases only if species not already found
                    if orf_match_level != "species":
                        if any(alias.lower() in line.lower() for alias in genus_aliases):
                            orf_match_level = "genus"
        except Exception as e:
            print(f"Warning: could not parse {orf_file}: {e}")

    return bin_match_level, bin_support, orf_match_level


def main(args):

    # Read expected species list
    expected_map = prep_expected_taxonomy(args.expected_csv)

    # build df from expected map
    expected_df = pd.DataFrame([
        {"match_name": acc,
        "expected_species": v["species"],
        "expected_genus": v["genus"]}
        for acc, v in expected_map.items()
    ])
    
    # Load bins file list
    bins_per_acc = load_bins_file(args.bins)

    # Read raw metagenome manysearch results
    results_df = pd.read_csv(args.mgx_manysearch_csv)

    # merge expected df with manysearch results df
    # this way we have all accessions, including any without manysearch matches
    results_df = expected_df.merge(
        pd.read_csv(args.mgx_manysearch_csv),
        on="match_name",
        how="left"
    )
   
    # Add columns for bin information
    results_df["total_n_bins"] = results_df["match_name"].map(lambda acc: bins_per_acc.get(acc, 0))
    results_df["no_bins"] = results_df["total_n_bins"] == 0

    # Add expected species and genus columns
    results_df["expected_species"] = None
    results_df["mgx_expected_species_found"] = False
    results_df["mgx_expected_genus_found"] = False

     # === Initialize bin columns ===

    # search bins x exact query genomes
    results_df["exact_bin_match_level"] = "unmatched"
    results_df["bin_ani_to_top_match_exact"] = None
    results_df["bin_ani_top_match_name"] = None
    results_df["exact_n_bins_with_match"] = 0

    # search bins x NCBI database
    results_df["ncbi_bin_match_level"] = "unmatched"
    results_df["bin_ani_to_top_match_ncbi"] = None
    results_df["bin_ani_top_match_name_ncbi"] = None
    results_df["ncbi_n_bins_with_match"] = 0

    # BAT results
    results_df["bat_bin_match_level"] = "unmatched"
    results_df["bat_bin_support"] = None
    results_df["bat_orf_match_level"] = "unmatched"


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
                results_df.at[idx, "exact_bin_match_level"] = "NA"
                results_df.at[idx, "ncbi_bin_match_level"] = "NA"
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
                results_df.at[idx, "exact_bin_match_level"] = "species"
                results_df.at[idx, "exact_n_bins_with_match"] = len(species_matches)
                results_df.at[idx, "bin_ani_to_top_match_exact"] = species_matches["average_containment_ani"].max()

            elif not genus_matches.empty:
                results_df.at[idx, "exact_bin_match_level"] = "genus"
                results_df.at[idx, "exact_n_bins_with_match"] = len(genus_matches)
                # If species match is False but genus match exists, best ANI is from genus matches
                results_df.at[idx, "bin_ani_to_top_match_exact"] = genus_matches["average_containment_ani"].max()
            
            # don't check anything we don't have in the database (e.g., viruses)
            if exp_genus and exp_genus.strip().lower() in NCBI_EXCLUDE_GENUS:
                results_df.at[idx, "ncbi_bin_match_level"] = "unmatched"
                # leave counts/ANI at defaults
                continue
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
                    # results_df.at[idx, "ncbi_bin_species_match"] = True
                    results_df.at[idx, "ncbi_bin_match_level"] = "species"
                    results_df.at[idx, "ncbi_n_bins_with_match"] = len(species_rows_ncbi)
                    results_df.at[idx, "bin_ani_to_top_match_ncbi"], results_df.at[idx, "bin_ani_top_match_name_ncbi"] = get_top_match_info(species_rows_ncbi)

                elif not genus_rows_ncbi.empty:
                    # results_df.at[idx, "ncbi_bin_genus_match"] = True
                    results_df.at[idx, "ncbi_bin_match_level"] = "genus"
                    results_df.at[idx, "ncbi_n_bins_with_match"] = len(genus_rows_ncbi)
                    results_df.at[idx, "bin_ani_to_top_match_ncbi"], results_df.at[idx, "bin_ani_top_match_name_ncbi"] = get_top_match_info(genus_rows_ncbi)

            # --- Check BAT results ---
            bat_bin_level, bat_bin_support, bat_orf_level = check_bat_for_accession(
                accession, exp_genus, exp_species, bat_dir="output.BAT"
            )
            results_df.at[idx, "bat_bin_match_level"] = bat_bin_level
            results_df.at[idx, "bat_bin_support"] = bat_bin_support
            results_df.at[idx, "bat_orf_match_level"] = bat_orf_level


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
        print("\nAll raw metagenome accessions had genome search matches.")

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
    bin_species_count = (valid_bin_rows["exact_bin_match_level"] == "species").sum()
    bin_genus_count = valid_bin_rows["exact_bin_match_level"].isin(["genus", "species"]).sum()

    print(f"Metagenomes with bin species match: {bin_species_count} / {len(valid_bin_rows)}")
    print(f"Metagenomes with bin genus match: {bin_genus_count} / {len(valid_bin_rows)}")

    # === Report missing bin matches ===
    missing_bin_matches = valid_bin_rows[valid_bin_rows["exact_bin_match_level"] == "unmatched"]

    if not missing_bin_matches.empty:
        print("No bin species/genus match:")
        for _, row in missing_bin_matches.iterrows():
            acc = row["match_name"]
            n_bins = row["total_n_bins"]
            print(f"  {acc} — expected {expected_map[acc]['species']} ({n_bins} bins)")
            bin_info = bins_df[bins_df["metagenome_accession"] == acc]
            if not bin_info.empty:
                print(bin_info.to_string(index=False))

    print("\n=== BINS X NCBI GENOMES ===")
    ncbi_species_count = (valid_bin_rows["ncbi_bin_match_level"] == "species").sum()
    ncbi_genus_count = valid_bin_rows["ncbi_bin_match_level"].isin(["genus", "species"]).sum()
    print(f"Metagenomes with bin species match: {ncbi_species_count} / {len(valid_bin_rows)}")
    print(f"Metagenomes with bin genus match: {ncbi_genus_count} / {len(valid_bin_rows)}")

    missing_ncbi = valid_bin_rows[valid_bin_rows["ncbi_bin_match_level"] == "unmatched"]

    if not missing_ncbi.empty:
        print("\nNo bin species/genus match (NCBI manysearch):")
        for _, r in missing_ncbi.iterrows():
            print(f"  {r['match_name']} — expected {expected_map[r['match_name']]['species']} ({r['total_n_bins']} bins)")

    print("\n=== BAT CLASSIFICATION SUMMARY ===")
    bat_species_count = (results_df["bat_bin_match_level"] == "species").sum()
    bat_genus_count = results_df["bat_bin_match_level"].isin(["genus", "species"]).sum()
    print(f"Metagenomes with BAT bin species match: {bat_species_count} / {len(results_df)}")
    print(f"Metagenomes with BAT bin genus match: {bat_genus_count} / {len(results_df)}")

    orf_species_count = (results_df["bat_orf_match_level"] == "species").sum()
    orf_genus_count = results_df["bat_orf_match_level"].isin(["genus", "species"]).sum()
    print(f"Metagenomes with BAT ORF species match: {orf_species_count} / {len(results_df)}")
    print(f"Metagenomes with BAT ORF genus match: {orf_genus_count} / {len(results_df)}")

    # --- Report missing BAT bin matches ---
    missing_bat_bins = results_df[results_df["bat_bin_match_level"] == "unmatched"]
    if not missing_bat_bins.empty:
        print("\nNo BAT bin species/genus match:")
        for _, row in missing_bat_bins.iterrows():
            acc = row["match_name"]
            n_bins = row["total_n_bins"]
            exp_species = expected_map.get(acc, {}).get("species", "N/A")
            print(f"  {acc} — expected {exp_species} ({n_bins} bins)")

    # --- Report missing BAT ORF matches ---
    missing_bat_orfs = results_df[results_df["bat_orf_match_level"] == "unmatched"]
    if not missing_bat_orfs.empty:
        print("\nNo BAT ORF species/genus match:")
        for _, row in missing_bat_orfs.iterrows():
            acc = row["match_name"]
            exp_species = expected_map.get(acc, {}).get("species", "N/A")
            print(f"  {acc} — expected {exp_species}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check metagenome search results against expected species list.")
    parser.add_argument("--mgx-manysearch-csv", help="CSV file with manysearch results.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.csv")
    parser.add_argument("--bin-manysearch-ncbi-csv", help="CSV file with bin search results.", default="output.manysearch/bins-x-ncbi-entire.manysearch.csv")
    parser.add_argument("--expected-csv", help="CSV file with expected metagenome accessions and species", default="multi/multi_mapping.csv")
    parser.add_argument("-o", "--output", help="Output CSV with added expected species/genus columns.", default="multi-aggregated-results.csv")
    parser.add_argument("--bins", help="path to all bins", default="multi_bins.txt")
    parser.add_argument("--multisearch-bins", help = "CSV file with multisearch bins", default="output.manysearch/bins-x-search-genomes.multisearch.default-thresh.csv")
    args = parser.parse_args()

    main(args)
