import argparse
import os, re
import csv
from dataclasses import dataclass, asdict
from typing import Optional
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

@dataclass
class AccessionResult:

    accession: str
    expected_genus: str
    expected_species: str
    
    # NCBI metadata
    assay_type: Optional[str] = None
    bioproject: Optional[str] = None
    biosample_link: Optional[str] = None
    collection_date_sam: Optional[str] = None
    geo_loc_name_country_calc: Optional[str] = None
    lat_lon: Optional[str] = None
    organism: Optional[str] = None

    # mgx results
    mgx_expected_genus_found: bool = False
    mgx_expected_species_found: bool = False
    mgx_containment: Optional[float] = None
    mgx_containment_target_in_query: Optional[float] = None
    mgx_f_weighted_target_in_query: Optional[float] = None
    
    # branchwater-web results
    mgx_k21_cANI: Optional[float] = None
    mgx_k21_containment: Optional[float] = None

    # bin results
    total_n_bins: int = 0
    exact_bin_match_level: str = "unmatched"
    exact_n_bins_with_match: int = 0
    bin_ani_to_top_match_exact: Optional[float] = None

    ncbi_bin_match_level: str = "unmatched"
    ncbi_n_bins_with_match: int = 0
    bin_ani_to_top_match_ncbi: Optional[float] = None
    bin_ani_top_match_name_ncbi: Optional[str] = None

    # BAT results
    bat_bin_match_level: str = "unmatched"
    bat_bin_support: Optional[float] = None
    bat_orf_match_level: str = "unmatched"

    # sendsketch results
    sendsketch_bin_match_level: str = "unmatched"
    

    def to_dict(self):
        return asdict(self)

def write_results_csv(results: list[AccessionResult], outpath: str):
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(AccessionResult.__annotations__.keys()))
        writer.writeheader()
        for res in results:
            writer.writerow(res.to_dict())

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


def summarize_results(results: list[AccessionResult], expected_map: dict):
    total_expected = len(expected_map)

    # === MGX ===
    found_accessions = {r.accession for r in results if r.mgx_expected_genus_found or r.mgx_expected_species_found}
    missing = [(acc, expected_map[acc]["species"]) for acc in expected_map if acc not in found_accessions]
    if missing:
        print("\nMissing accessions:")
        for acc, sp in missing:
            print(f"  {acc} — expected {sp}")
    else:
        print("\nAll raw metagenome accessions had genome search matches.")

    correct_species = sum(r.mgx_expected_species_found for r in results)
    correct_genus = sum(r.mgx_expected_genus_found for r in results)
    print("\n=== SEARCH GENOMES X METAGENOMES ===")
    print(f"Total expected metagenome matches: {total_expected}")
    print(f"Correct species matches: {correct_species} / {total_expected}")
    print(f"Correct genus matches: {correct_genus} / {total_expected}")

    # === Exact bins ===
    valid_bins = [r for r in results if r.total_n_bins > 0]
    bin_species = sum(r.exact_bin_match_level == "species" for r in valid_bins)
    bin_genus = sum(r.exact_bin_match_level in ("species","genus") for r in valid_bins)
    print("\n=== BINS X SEARCH GENOMES ===")
    print(f"Metagenomes with bin species match: {bin_species} / {len(valid_bins)}")
    print(f"Metagenomes with bin genus match: {bin_genus} / {len(valid_bins)}")

    missing_bin_matches = [r for r in valid_bins if r.exact_bin_match_level == "unmatched"]
    if missing_bin_matches:
        print("No bin species/genus match:")
        for r in missing_bin_matches:
            print(f"  {r.accession} — expected {r.expected_species} ({r.total_n_bins} bins)")

    # === NCBI bins ===
    ncbi_species = sum(r.ncbi_bin_match_level == "species" for r in valid_bins)
    ncbi_genus = sum(r.ncbi_bin_match_level in ("species","genus") for r in valid_bins)
    print("\n=== BINS X NCBI DATABASE ===")
    print(f"Metagenomes with bin species match: {ncbi_species} / {len(valid_bins)}")
    print(f"Metagenomes with bin genus match: {ncbi_genus} / {len(valid_bins)}")

    missing_ncbi = [r for r in valid_bins if r.ncbi_bin_match_level == "unmatched"]
    if missing_ncbi:
        print("\nNo bin species/genus match (NCBI manysearch):")
        for r in missing_ncbi:
            print(f"  {r.accession} — expected {r.expected_species} ({r.total_n_bins} bins)")

    # === BAT ===
    bat_species = sum(r.bat_bin_match_level == "species" for r in results)
    bat_genus = sum(r.bat_bin_match_level in ("species","genus") for r in results)
    print("\n=== BAT CLASSIFICATION SUMMARY ===")
    print(f"Metagenomes with BAT bin species match: {bat_species} / {len(results)}")
    print(f"Metagenomes with BAT bin genus match: {bat_genus} / {len(results)}")

    orf_species = sum(r.bat_orf_match_level == "species" for r in results)
    orf_genus = sum(r.bat_orf_match_level in ("species","genus") for r in results)
    print(f"Metagenomes with BAT ORF species match: {orf_species} / {len(results)}")
    print(f"Metagenomes with BAT ORF genus match: {orf_genus} / {len(results)}")

    missing_bat_bins = [r for r in results if r.bat_bin_match_level == "unmatched"]
    if missing_bat_bins:
        print("\nNo BAT bin species/genus match:")
        for r in missing_bat_bins:
            print(f"  {r.accession} — expected {r.expected_species} ({r.total_n_bins} bins)")

    missing_bat_orfs = [r for r in results if r.bat_orf_match_level == "unmatched"]
    if missing_bat_orfs:
        print("\nNo BAT ORF species/genus match:")
        for r in missing_bat_orfs:
            print(f"  {r.accession} — expected {r.expected_species}")
    
    # === SENDSKETCH ===
    valid_ss = [r for r in results if r.sendsketch_bin_match_level != "NA"]

    ss_species = sum(r.sendsketch_bin_match_level == "species" for r in valid_ss)
    ss_genus   = sum(r.sendsketch_bin_match_level in ("species", "genus") for r in valid_ss)

    print("\n=== SENDSKETCH BINS ===")
    print(f"Metagenomes with SendSketch species match: {ss_species} / {len(valid_ss)}")
    print(f"Metagenomes with SendSketch genus match: {ss_genus} / {len(valid_ss)}")

    missing_ss = [r for r in valid_ss if r.sendsketch_bin_match_level == "unmatched"]
    if missing_ss:
        print("\nNo SendSketch species/genus match:")
        for r in missing_ss:
            print(f"  {r.accession} — expected {r.expected_species}")


def accession_summary_table(results: list[AccessionResult]) -> list[dict]:
    """
    Build a table of accession-level summary results across all tools.
    Returns list of dicts (rows).
    """
    table = []
    for r in results:
        # MGX: grab k21 cANI, % metagenome
        mgx_cANI = (
            f"{float(r.mgx_k21_cANI) * 100:.1f}"
            if r.mgx_k21_cANI not in (None, "") else ""
        )
        mgx_pct = (
            f"{float(r.mgx_f_weighted_target_in_query) * 100:.1f}"
            if r.mgx_f_weighted_target_in_query not in (None, "") else ""
        )

        bat_support_pct = (
            f"{float(r.bat_bin_support) * 100:.1f}"
            if r.bat_bin_support not in (None, "") else ""
        )


        row = {
            "Accession": r.accession,
            "ExpectedSpecies": r.expected_species,
            "cANI": mgx_cANI,
            "% mgx": mgx_pct,
            "Bins_x_SearchGx": r.exact_bin_match_level,
            "Bins_x_NCBI": r.ncbi_bin_match_level,
            "BAT_Bins": r.bat_bin_match_level,
            "BAT_Bins_support": bat_support_pct,
            "BAT_ORFs": r.bat_orf_match_level,
            "SendSketch": r.sendsketch_bin_match_level,
        }
        table.append(row)
    return table

def write_accession_summary(results: list[AccessionResult],
                            out_csv: str | None = None):
    """
    Write metagenome-level summary to stdout, and also to CSV if out_csv is provided.
    """
    table = accession_summary_table(results)
    header = list(table[0].keys()) if table else []

    # --- Print to stdout ---
    print("\n=== METAGENOME-LEVEL SUMMARY ===")
    print("\t".join(header))
    for row in table:
        print("\t".join(str(row[h]) for h in header))

    # --- Write CSV if requested ---
    if out_csv:
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(table)
        print(f"\nWrote accession summary to '{out_csv}'")


def prep_expected_taxonomy(expected_csv: str) -> dict:
    """
    Load expected taxonomy from CSV and prepare a mapping for quick lookup.
    CSV format: accession, expected_genus, expected_species
    Returns OrderedDict like:
      {
        "SRR12345": {"genus": "Aspergillus", "species": "Aspergillus sydowii"},
        ...
      }
    """
    expected_map = OrderedDict()
    with open(expected_csv, encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or len(row) < 3:
                continue
            accession, genus, species = row[0].strip(), row[1].strip(), row[2].strip()
            accession = accession.lstrip("\ufeff")
            species = species.replace("_", " ")  # underscores → spaces
            expected_map[accession] = {
                "genus": genus,
                "species": f"{genus} {species}"
            }
    return expected_map

def load_bins_file(bins_file: str) -> dict:
    """
    Load bins file and return a mapping of accession to number of bins.
    """
    bins_per_acc = {}
    with open(bins_file, encoding="utf-8-sig") as f:
        for line in f:
            path = line.strip()
            if not path:
                continue
            acc = path.split("/")[2]  # e.g., SRR11734785 from multi/bins/SRR11734785/file.fa
            bins_per_acc[acc] = bins_per_acc.get(acc, 0) + 1
    return bins_per_acc

def load_mgx_results(path: str) -> dict:
    """
    Return accession -> row (full dict) from mgx manysearch CSV.
    Uses 'match_name' as accession key.
    """
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["match_name"].strip()
            out[acc] = row
    return out

def load_sendsketch_results(path: str) -> dict:
    """
    Load SendSketch results/metadata CSV into accession -> row dict.
    Keys are 'acc' (metagenome accession).
    """
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["acc"].strip()
            out[acc] = row
    return out


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


def check_exact_bins_for_accession(acc: str,
                                   exp_genus: str,
                                   exp_species: str,
                                   bins_df: list[dict]) -> tuple[str, int, float | None]:
    """
    Return (match_level, n_bins_with_match, best_ani).
    bins_df should be a list of dict rows from the bin multisearch CSV,
    with at least: 'metagenome_accession', 'match_name', 'average_containment_ani'.
    """
    g_aliases, s_aliases = get_aliases(exp_genus, exp_species)

    rows = [r for r in bins_df if r["metagenome_accession"] == acc]
    if not rows:
        return "unmatched", 0, None

    # species matches
    s_rows = [r for r in rows if any(a.lower() in str(r["match_name"]).lower() for a in s_aliases)]
    g_rows = [r for r in rows if any(a.lower() in str(r["match_name"]).lower() for a in g_aliases)]

    if s_rows:
        best = max(r.get("average_containment_ani", 0) for r in s_rows)
        return "species", len(s_rows), best
    elif g_rows:
        best = max(r.get("average_containment_ani", 0) for r in g_rows)
        return "genus", len(g_rows), best
    else:
        return "unmatched", 0, None


def check_ncbi_bins_for_accession(acc: str,
                                  exp_genus: str,
                                  exp_species: str,
                                  ncbi_bins: list[dict]) -> tuple[str, int, float | None, str | None]:
    """
    Return (match_level, n_bins_with_match, best_ani, best_match_name).
    ncbi_bins should be a list of dict rows from bin-manysearch-ncbi CSV.
    """
    g_aliases, s_aliases = get_aliases(exp_genus, exp_species)

    rows = [r for r in ncbi_bins if r["metagenome_accession"] == acc]
    if not rows:
        return "unmatched", 0, None, None

    # species first
    s_rows = [r for r in rows if any(a.lower() in str(r["match_name"]).lower() for a in s_aliases)]
    g_rows = [r for r in rows if any(a.lower() in str(r["match_name"]).lower() for a in g_aliases)]

    if s_rows:
        best_row = max(s_rows, key=lambda r: r.get("average_containment_ani", 0))
        return "species", len(s_rows), best_row.get("average_containment_ani"), best_row.get("match_name")
    elif g_rows:
        best_row = max(g_rows, key=lambda r: r.get("average_containment_ani", 0))
        return "genus", len(g_rows), best_row.get("average_containment_ani"), best_row.get("match_name")
    else:
        return "unmatched", 0, None, None

def process_bat_bin_records(records: list[dict], g_aliases: list[str], s_aliases: list[str]) -> tuple[str, Optional[float]]:
    """
    Process parsed BAT bin2classification.taxnames.txt rows.
    Return (bin_match_level, bin_support).
    """
    bin_match_level = "unmatched"
    bin_support = None

    for row in records:
        lineage_text = row.get("lineage scores (f: 0.30)", "")  # contains "Genus: 0.95 Species: 0.91" style
        genus = row['genus']
        species = row['species']

        if species != "no support" and species != "NA" and species != None:
            try:
                name, score_str = species.split(": ")
                # account for any trailing '*'
                name = name.rstrip("*").strip()
                score = float(score_str)
            except ValueError:
                print(f"WARNING: could not parse species and score from row {row}")
                continue
            if any(alias.lower() == name.lower() for alias in s_aliases):
                # make sure we're saving the best species score
                if bin_match_level != "species" or (bin_support is None or score > bin_support):
                    bin_match_level, bin_support = "species", score

        # --- Genus check (only if no species support was found) ---
        if bin_match_level != "species":
            # check genus aliases
            if genus != "no support" and genus != "NA" and species != None:
                try:
                    name, score_str = genus.split(": ")
                    # account for any trailing '*'
                    name = name.rstrip("*").strip()
                    score = float(score_str)
                except ValueError:
                    print(f"WARNING: could not parse genus and score from row {row}")
                    continue
                if any(alias.lower() == name.lower() for alias in g_aliases):
                    if bin_match_level != "genus" or (bin_support is None or score > bin_support):
                        bin_match_level, bin_support = "genus", score

    return bin_match_level, bin_support

def load_bat_orf_file(path: str) -> list[dict]:
    """
    Load a BAT ORF2LCA.taxnames.txt file into dicts.
    First 5 fields are fixed, the rest are merged into 'full_lineage_names'.
    """
    records = []
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        # Normalize to expected first 5 columns
        base_fields = ["ORF", "bin", "n_hits", "lineage", "top_bit_score"]
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            base = row[:5]
            extra = row[5:]
            rec = {k: v for k, v in zip(base_fields, base)}
            rec["full_lineage_names"] = "\t".join(extra)
            records.append(rec)
    return records


def process_bat_orf_records(records: list[dict], g_aliases: list[str], s_aliases: list[str]) -> str:
    """
    Process parsed BAT ORF2LCA.taxnames.txt rows.
    Return orf_match_level ("species", "genus", "unmatched").
    """
    orf_match_level = "unmatched"

    for row in records:
        lineage_text = row.get("full_lineage_names", "")
        if not lineage_text:
            continue

        if any(alias.lower() in lineage_text.lower() for alias in s_aliases):
            return "species"
        if orf_match_level != "species" and any(alias.lower() in lineage_text.lower() for alias in g_aliases):
            orf_match_level = "genus"

    return orf_match_level


def check_bat_for_accession(acc: str,
                            exp_genus: str,
                            exp_species: str,
                            bat_dir="output.BAT") -> tuple[str, Optional[float], str]:
    """
    Return (bat_bin_match_level, bat_bin_support, bat_orf_match_level).
    """
    acc_dir = os.path.join(bat_dir, acc)
    bin_file = os.path.join(acc_dir, "out.BAT.bin2classification.taxnames.txt")
    orf_file = os.path.join(acc_dir, "out.BAT.ORF2LCA.taxnames.txt")

    g_aliases, s_aliases = get_aliases(exp_genus, exp_species)

    # --- Bin classifications ---
    if os.path.exists(bin_file):
        with open(bin_file, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            bin_records = list(reader)
        bin_match_level, bin_support = process_bat_bin_records(bin_records, g_aliases, s_aliases)
    else:
        bin_match_level, bin_support = "NA", None

    # --- ORF classifications ---
    if os.path.exists(orf_file):
        recs = load_bat_orf_file(orf_file)
        orf_match_level = process_bat_orf_records(recs, g_aliases, s_aliases)
    else:
        orf_match_level = "NA"

    return bin_match_level, bin_support, orf_match_level


def main(args):
    # --- Load expected taxonomy ---
    expected_map = prep_expected_taxonomy(args.expected_csv)

    # --- Load bins count ---
    bins_per_acc = load_bins_file(args.bins)

    # --- Load MGX results ---
    mgx_results = load_mgx_results(args.mgx_manysearch_csv)

    # --- Load SendSketch results ---
    sendsketch_results = load_sendsketch_results(args.sendsketch_csv)

    # --- Load bin search results ---
    with open(args.multisearch_bins, newline="") as f:
        bins_df = list(csv.DictReader(f))
    
    # derive metagenome accession (SRRxxxx from SRRxxxx.y)
    for r in bins_df:
        if "query_name" in r and r["query_name"]:
            r["metagenome_accession"] = r["query_name"].split(".", 1)[0]

    with open(args.bin_manysearch_ncbi_csv, newline="") as f:
        ncbi_bins = list(csv.DictReader(f))
    
    # derive metagenome accession (SRRxxxx from SRRxxxx.y)
    for r in ncbi_bins:
        if "query_name" in r and r["query_name"]:
            r["metagenome_accession"] = r["query_name"].split(".", 1)[0]

    # --- Process each accession ---
    results: list[AccessionResult] = []
    for acc, tax in expected_map.items():
        res = AccessionResult(accession=acc,
                              expected_genus=tax["genus"],
                              expected_species=tax["species"])

        res.total_n_bins = bins_per_acc.get(acc, 0)

        # MGX Results
        mgx_row = mgx_results.get(acc)
        if mgx_row:
            res.mgx_expected_genus_found, res.mgx_expected_species_found = \
                check_expected_in_text(mgx_row["query_name"], res.expected_genus, res.expected_species)
            res.mgx_containment = float(mgx_row.get("containment") or 0)
            res.mgx_containment_target_in_query = float(mgx_row.get("containment_target_in_query") or 0)
            res.mgx_f_weighted_target_in_query = float(mgx_row.get("f_weighted_target_in_query") or 0)
        
        # Bin Results
        
        # Sendsketch + metadata
        sendsketch_row = sendsketch_results.get(acc)
        if sendsketch_row:
            # first, get metadata
            res.assay_type = sendsketch_row.get("assay_type")
            res.bioproject = sendsketch_row.get("bioproject")
            res.biosample_link = sendsketch_row.get("biosample_link")
            res.collection_date_sam = sendsketch_row.get("collection_date_sam")
            res.geo_loc_name_country_calc = sendsketch_row.get("geo_loc_name_country_calc")
            res.lat_lon = sendsketch_row.get("lat_lon")
            res.organism = sendsketch_row.get("organism")

            # now branchwater web mgx search results
            res.mgx_k21_containment = sendsketch_row.get("containment")
            res.mgx_k21_cANI = sendsketch_row.get("cANI")

            # now actual sendsketch aggregated results
            res.sendsketch_bin_match_level = sendsketch_row.get("match level", "unmatched")


        # sourmash x exact search genomes
        res.exact_bin_match_level, res.exact_n_bins_with_match, res.bin_ani_to_top_match_exact = \
            check_exact_bins_for_accession(acc, res.expected_genus, res.expected_species, bins_df)

        # sourmash search x NCBI
        res.ncbi_bin_match_level, res.ncbi_n_bins_with_match, \
        res.bin_ani_to_top_match_ncbi, res.bin_ani_top_match_name_ncbi = \
            check_ncbi_bins_for_accession(acc, res.expected_genus, res.expected_species, ncbi_bins)

        # BAT annotation
        res.bat_bin_match_level, res.bat_bin_support, res.bat_orf_match_level = \
            check_bat_for_accession(acc, res.expected_genus, res.expected_species, "output.BAT")
        
        # --- enforce NA if no bins ---
        if res.total_n_bins == 0:
            res.exact_bin_match_level = "NA"
            res.ncbi_bin_match_level = "NA"
            res.bat_bin_match_level = "NA"
            res.bat_orf_match_level = "NA"

        results.append(res)

        # fix virus genus name
        if res.expected_genus == "Acanthamoeba":
            res.expected_genus = "Mimivirus"

    # --- Summarize by tool ---
    summarize_results(results, expected_map)
    
    # --- Summarize by Accession ---
    write_accession_summary(results, out_csv="multi-accession-summary.csv")
    
    # --- Write Full Aggregated CSV ---
    print(f"Writing full results to '{args.output}'...")
    write_results_csv(results, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate bin search and annotation results.")
    parser.add_argument("--expected-csv", help="CSV file with expected metagenome accessions and species", default="multi/multi_mapping.csv")
    parser.add_argument("--mgx-manysearch-csv", help="CSV file with manysearch results.", default="output.manysearch/search-genomes-x-brmetagenomes.manysearch.csv")
    parser.add_argument("--bin-manysearch-ncbi-csv", help="CSV file with bin search results.", default="output.manysearch/bins-x-ncbi-entire.manysearch.csv")
    parser.add_argument("-o", "--output", help="Output CSV with added expected species/genus columns.", default="multi-aggregated-results.csv")
    parser.add_argument("--bins", help="path to all bins", default="multi_bins.txt")
    parser.add_argument("--multisearch-bins", help = "CSV file with multisearch bins", default="output.manysearch/bins-x-search-genomes.multisearch.default-thresh.csv")
    parser.add_argument("--sendsketch-csv", help="CSV file with SendSketch results.", default="multi/multi.sendsketch.csv")

    args = parser.parse_args()

    main(args)
