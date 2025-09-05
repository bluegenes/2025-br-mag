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
    # original branchwater-web results
    mgx_k21_cANI: Optional[float] = None
    mgx_k21_containment: Optional[float] = None
    # k31 manysearch results
    mgx_expected_genus_found: bool = False
    mgx_expected_species_found: bool = False
    mgx_k31_containment: Optional[float] = None
    mgx_k31_containment_target_in_query: Optional[float] = None
    mgx_k31_f_weighted_target_in_query: Optional[float] = None

    # bin results
    nbins: int = 0
    bins_x_searchgx_match_level: str = "unmatched"
    bins_x_searchgx_nbins_matched: int = 0
    bins_x_searchgx_topANI: Optional[float] = None

    bins_x_ncbi_match_level: str = "unmatched"
    bins_x_ncbi_nbins_matched: int = 0
    bins_x_ncbi_topANI: Optional[float] = None
    bins_x_ncbi_top_match_name: Optional[str] = None

    # BAT results
    bins_x_batNR_match_level: str = "unmatched"
    bins_x_batNR_support: Optional[float] = None
    orfs_x_batNR_match_level: str = "unmatched"

    # sendsketch results
    bins_x_sendsketch_match_level: str = "unmatched"


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
    valid_bins = [r for r in results if r.nbins > 0]
    bin_species = sum(r.bins_x_searchgx_match_level == "species" for r in valid_bins)
    bin_genus = sum(r.bins_x_searchgx_match_level in ("species","genus") for r in valid_bins)
    print("\n=== BINS X SEARCH GENOMES ===")
    print(f"Metagenomes with bin species match: {bin_species} / {len(valid_bins)}")
    print(f"Metagenomes with bin genus match: {bin_genus} / {len(valid_bins)}")

    missing_bin_matches = [r for r in valid_bins if r.bins_x_searchgx_match_level == "unmatched"]
    if missing_bin_matches:
        print("No bin species/genus match:")
        for r in missing_bin_matches:
            print(f"  {r.accession} — expected {r.expected_species} ({r.nbins} bins)")

    # === NCBI bins ===
    ncbi_species = sum(r.bins_x_ncbi_match_level == "species" for r in valid_bins)
    ncbi_genus = sum(r.bins_x_ncbi_match_level in ("species","genus") for r in valid_bins)
    print("\n=== BINS X NCBI DATABASE ===")
    print(f"Metagenomes with bin species match: {ncbi_species} / {len(valid_bins)}")
    print(f"Metagenomes with bin genus match: {ncbi_genus} / {len(valid_bins)}")

    missing_ncbi = [r for r in valid_bins if r.bins_x_ncbi_match_level == "unmatched"]
    if missing_ncbi:
        print("\nNo bin species/genus match (NCBI manysearch):")
        for r in missing_ncbi:
            print(f"  {r.accession} — expected {r.expected_species} ({r.nbins} bins)")

    # === BAT ===
    bat_species = sum(r.bins_x_batNR_match_level == "species" for r in results)
    bat_genus = sum(r.bins_x_batNR_match_level in ("species","genus") for r in results)
    print("\n=== BAT CLASSIFICATION SUMMARY ===")
    print(f"Metagenomes with BAT bin species match: {bat_species} / {len(results)}")
    print(f"Metagenomes with BAT bin genus match: {bat_genus} / {len(results)}")

    orf_species = sum(r.orfs_x_batNR_match_level == "species" for r in results)
    orf_genus = sum(r.orfs_x_batNR_match_level in ("species","genus") for r in results)
    print(f"Metagenomes with BAT ORF species match: {orf_species} / {len(results)}")
    print(f"Metagenomes with BAT ORF genus match: {orf_genus} / {len(results)}")

    missing_bat_bins = [r for r in results if r.bins_x_batNR_match_level == "unmatched"]
    if missing_bat_bins:
        print("\nNo BAT bin species/genus match:")
        for r in missing_bat_bins:
            print(f"  {r.accession} — expected {r.expected_species} ({r.nbins} bins)")

    missing_bat_orfs = [r for r in results if r.orfs_x_batNR_match_level == "unmatched"]
    if missing_bat_orfs:
        print("\nNo BAT ORF species/genus match:")
        for r in missing_bat_orfs:
            print(f"  {r.accession} — expected {r.expected_species}")

    # === SENDSKETCH ===
    valid_ss = [r for r in results if r.bins_x_sendsketch_match_level != "NA"]

    ss_species = sum(r.bins_x_sendsketch_match_level == "species" for r in valid_ss)
    ss_genus   = sum(r.bins_x_sendsketch_match_level in ("species", "genus") for r in valid_ss)

    print("\n=== SENDSKETCH BINS ===")
    print(f"Metagenomes with SendSketch species match: {ss_species} / {len(valid_ss)}")
    print(f"Metagenomes with SendSketch genus match: {ss_genus} / {len(valid_ss)}")

    missing_ss = [r for r in valid_ss if r.bins_x_sendsketch_match_level == "unmatched"]
    if missing_ss:
        print("\nNo SendSketch species/genus match:")
        for r in missing_ss:
            print(f"  {r.accession} — expected {r.expected_species}")


def accession_summary_table(results: list[AccessionResult]) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    """
    Build a table of accession-level summary results across all tools.
    Returns (all_rows, unmatched_rows).
    """
    table = []
    nobins = []
    unmatched = []
    orf_match_only = []
    for r in results:
        mgx_cANI = (
            f"{float(r.mgx_k21_cANI) * 100:.1f}"
            if r.mgx_k21_cANI not in (None, "") else ""
        )
        mgx_pct = (
            f"{float(r.mgx_k31_f_weighted_target_in_query) * 100:.1f}"
            if r.mgx_k31_f_weighted_target_in_query not in (None, "") else ""
        )
        bat_support_pct = (
            f"{float(r.bins_x_batNR_support) * 100:.1f}"
            if r.bins_x_batNR_support not in (None, "") else ""
        )

        row = {
            "Accession": r.accession,
            "ExpectedSpecies": r.expected_species,
            "cANI": mgx_cANI,
            "% mgx": mgx_pct,
            "n_bins": r.nbins,
            "Bins_x_SearchGx": r.bins_x_searchgx_match_level,
            "Bins_x_NCBI": r.bins_x_ncbi_match_level,
            "BAT_Bins": r.bins_x_batNR_match_level,
            "BAT_Bins_support": bat_support_pct,
            "BAT_ORFs": r.orfs_x_batNR_match_level,
            "SendSketch": r.bins_x_sendsketch_match_level,
        }

        # separate no bins
        if r.nbins == 0:
            nobins.append(row)
            continue

        # separate unmatched bins (no species/genus match in any tool)

        # check if unmatched across all tools
        bin_has_species = (
            r.bins_x_searchgx_match_level == "species"
            or r.bins_x_ncbi_match_level == "species"
            or r.bins_x_batNR_match_level == "species"
            or r.bins_x_sendsketch_match_level == "species"
        )
        bin_has_genus = (
            r.bins_x_searchgx_match_level in ("species", "genus")
            or r.bins_x_ncbi_match_level in ("species", "genus")
            or r.bins_x_batNR_match_level in ("species", "genus")
            or r.bins_x_sendsketch_match_level in ("species", "genus")
        )

        orf_has_match = r.orfs_x_batNR_match_level in ("species", "genus")

        if not bin_has_species and not bin_has_genus:
            if orf_has_match:
                orf_match_only.append(row)
                continue
            else:
                unmatched.append(row)
                continue

        table.append(row)

    return table, nobins, unmatched, orf_match_only

def write_accession_summary(results: list[AccessionResult],
                            out_csv: str | None = None):
    """
    Write metagenome-level summary to stdout, and also to CSV if out_csv is provided.
    """
    table, nobins, unmatched, orf_match_only = accession_summary_table(results)
    header = list(table[0].keys()) if table else []

    # --- Print matched summary ---
    print("\n=== BIN SUMMARY ===")
    print("\t".join(header))
    for row in table:
        print("\t".join(str(row[h]) for h in header))

    # --- Print no bins ---
    if nobins:
        print("\n=== NO BINS ===")
        print("\t".join(header))
        for row in nobins:
            print("\t".join(str(row[h]) for h in header))

    # --- Print ORF match only ---
    if orf_match_only:
        print("\n=== ORF MATCH ONLY (no bin-level species/genus match) ===")
        print("\t".join(header))
        for row in orf_match_only:
            print("\t".join(str(row[h]) for h in header))

    # --- Print unmatched ---
    if unmatched:
        print("\n=== UNMATCHED BINS (no species/genus match) ===")
        print("\t".join(header))
        for row in unmatched:
            print("\t".join(str(row[h]) for h in header))

    # --- Write CSV if requested ---
    if out_csv:
        with open(out_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(table)
            writer.writerows(nobins)
            writer.writerows(orf_match_only)
            writer.writerows(unmatched)
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

    # identify matches to the search genome
    match_rows = [r for r in rows if any(a.lower() in str(r["match_name"]).lower() for a in s_aliases)]

    # estimate species-level match: containment ANI in either direction is ANI >= 0.95; else genus
    if match_rows:
        valid_ani = [float(r.get("max_containment_ani", 0)) for r in match_rows if r.get("max_containment_ani") not in (None, "", "NA")]
        best = max(valid_ani) if valid_ani else 0
        if best >= 0.95:
            return "species", len(match_rows), best
        else:
            return "genus", len(match_rows), best
    else:
        return "unmatched", 0, None


def aggregate_sourmash_ncbi(acc: str,
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
        # special case: NCBI is missing Mimivirus/Acanthamoeba
        if exp_genus in ("Acanthamoeba", "Mimivirus") or "mimivirus" in exp_species.lower():
            return "NA", 0, None, None
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
        # again, apply special-case NA
        if exp_genus in ("Acanthamoeba", "Mimivirus") or "mimivirus" in exp_species.lower():
            return "NA", 0, None, None
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
    Return (bins_x_batNR_match_level, bins_x_batNR_support, orfs_x_batNR_match_level).
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

        res.nbins = bins_per_acc.get(acc, 0)

        # MGX Results
        mgx_row = mgx_results.get(acc)
        if mgx_row:
            res.mgx_expected_genus_found, res.mgx_expected_species_found = \
                check_expected_in_text(mgx_row["query_name"], res.expected_genus, res.expected_species)
            res.mgx_k31_containment = float(mgx_row.get("containment") or 0)
            res.mgx_k31_containment_target_in_query = float(mgx_row.get("containment_target_in_query") or 0)
            res.mgx_k31_f_weighted_target_in_query = float(mgx_row.get("f_weighted_target_in_query") or 0)

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
            res.bins_x_sendsketch_match_level = sendsketch_row.get("match level", "unmatched")
        # sourmash bins x exact search genomes
            res.bins_x_searchgx_match_level, res.bins_x_searchgx_nbins_matched, res.bins_x_searchgx_topANI = \
            check_exact_bins_for_accession(acc, res.expected_genus, res.expected_species, bins_df)

        # sourmash search bins x NCBI
        res.bins_x_ncbi_match_level, res.bins_x_ncbi_nbins_matched, \
        res.bins_x_ncbi_topANI, res.bins_x_ncbi_top_match_name = \
            aggregate_sourmash_ncbi(acc, res.expected_genus, res.expected_species, ncbi_bins)

        # BAT bin annotation
        res.bins_x_batNR_match_level, res.bins_x_batNR_support, res.orfs_x_batNR_match_level = \
            check_bat_for_accession(acc, res.expected_genus, res.expected_species, "output.BAT")

        # --- enforce NA if no bins ---
        if res.nbins == 0:
            res.bins_x_searchgx_match_level = "NA"
            res.bins_x_ncbi_match_level = "NA"
            res.bins_x_batNR_match_level = "NA"
            res.orfs_x_batNR_match_level = "NA"

        results.append(res)

        # fix virus genus name in output CSVs
        if res.expected_genus == "Acanthamoeba":
            res.expected_genus = "Mimivirus"

    # --- Summarize by tool ---
    summarize_results(results, expected_map)

    # --- Summarize by Accession ---
    write_accession_summary(results, out_csv=args.bin_summary)

    # --- Write Full Aggregated CSV ---
    print(f"Writing full results to '{args.output}'...")
    write_results_csv(results, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate bin search and annotation results.")
    parser.add_argument("--expected-csv", help="CSV file with expected metagenome accessions and species", default="multi/multi_mapping.csv")
    parser.add_argument("--mgx-manysearch-csv", help="CSV file with manysearch results.", default="output.sourmash/search-genomes-x-brmetagenomes.manysearch.csv")
    parser.add_argument("--bin-manysearch-ncbi-csv", help="CSV file with bin search results.", default="output.sourmash/bins-x-ncbi-euks.multisearch.sc1000.csv")
    parser.add_argument("-o", "--output", help="Output CSV with mgx, bin summarization.", default="multi.aggregated-results.csv")
    parser.add_argument("--bin-summary", help="Output CSV with only bin summary results.", default="multi.aggregated-results-bins.csv")
    parser.add_argument("--bins", help="path to all bins", default="multi_bins.txt")
    parser.add_argument("--multisearch-bins", help = "CSV file with multisearch bins", default="output.sourmash/bins-x-search-genomes.multisearch.csv")
    parser.add_argument("--sendsketch-csv", help="CSV file with SendSketch results.", default="multi/multispecies_tracking.csv")

    args = parser.parse_args()

    main(args)
