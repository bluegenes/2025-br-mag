import os
import pandas as pd
from pathlib import Path

# 1. classify bins using sendsketch, gtdbtk, and sourmash (...and anything else?)
    # --> gtdbtk doesn't do euks, find an alternate.
# 2. summarize the classification results
# 3. maybe run genome-grist to do mapping-based verification?

outdir = "output.classify_bins"
basename = "multi"


# Load bin information
bininfo = pd.read_csv("multi/bins.txt", header=None, names=["bin_file"])
# Extract sample/sample-bin names from the file names
# file name structure: MEGAHIT-MetaBAT2-{sample_bin}.fa.gz
SAMPLEBINS = bininfo.bin_file.str.extract(r'MEGAHIT-MetaBAT2-(?P<sample_bin>.+)\.fa\.gz')
SAMPLES = SAMPLEBINS.sample_bin.str.split('.', expand=True)[0].unique().tolist()
SAMPLE_BINS = SAMPLEBINS.sample_bin.tolist()


# some sample-bins don't exist??? check which
pattern = Path("multi/bins").glob("MEGAHIT-MetaBAT2-*.fa.gz")
existing_bins = [p.stem.replace("MEGAHIT-MetaBAT2-", "").replace(".fa", "")       # strip prefix
                    for p in pattern]

missing_bins = sorted(set(SAMPLE_BINS) - set(existing_bins))
print(f"{len(missing_bins)} bins are missing:\n")
for b in missing_bins:
    print(b)
# import pdb;pdb.set_trace()

# turn that into whatever containers you need
SAMPLE_BINS_EXISTING = existing_bins
SAMPLES_EXISTING      = {b.split('.')[0] for b in existing_bins}


# samples_txt = "multi/multinames.txt"
# SAMPLES = pd.read_csv(samples_txt, header=None)[0].tolist()
# multi_csv = pd.read_csv("multi/multi_mapping.csv", header=None)


rule all:
    input:
        expand(f"{outdir}/sourmash/{basename}-x-{{db}}.classifications.csv", db=["gtdb", "ncbi"]),
        #expand("results/sourmash/{sample}.taxonomy.csv", sample=SAMPLEBINS),
        # expand("results/sendsketch/{samplebin}.txt", sample=SAMPLEBINS),
        # expand("results/gtdbtk/{sample}/gtdbtk.bac120.summary.tsv", sample=SAMPLEBINS),
        # "results/summary/combined_taxonomy_summary.csv"

rule write_manysketch_file:
    input:
        fasta=expand("multi/bins/MEGAHIT-MetaBAT2-{sample_bin}.fa.gz", sample_bin=SAMPLE_BINS_EXISTING)
    output:
        csv=f"{outdir}/multi.manysketch.csv"
    run:
        # write file with name,genome_filename,protein_filename for all input files
        with open(output.csv, 'w') as f:
            f.write("name,genome_filename,protein_filename\n")
            # iterate over input fasta files and write to csv
            for fasta in input.fasta:
                sample_bin = os.path.basename(fasta).replace("MEGAHIT-MetaBAT2-", "").replace(".fa.gz", "")
                f.write(f"{sample_bin},{fasta},\n")

rule sourmash_sketch:
    input:
        csv=f"{outdir}/{basename}.manysketch.csv"
    output:
        sig=f"{outdir}/sourmash/{basename}.sig.zip"
    shell:
        """
        sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000 {input.csv} -o {output.sig}
        """

rule sourmash_fastmultigather_gtdb:
    input:
        zip=f"{outdir}/sourmash/{basename}.sig.zip",
        db= "/group/ctbrowngrp5/sourmash-db/gtdb-rs226/gtdb-rs226.k31.rocksdb"
    output:
        gather=f"{outdir}/sourmash/{basename}-x-gtdb.fmg.csv"
    shell:
        """
        sourmash scripts fastmultigather {input.zip} {input.db} -o {output.gather}
        """

rule sourmash_fastmultigather_ncbi:
    input:
        zip=f"{outdir}/sourmash/{basename}.sig.zip",
        db= "/group/ctbrowngrp5/sourmash-db/entire-2025-01-21/entire-2025-01-21.k51.rocksdb"
    output:
        gather=f"{outdir}/sourmash/{basename}-x-ncbi.fmg.csv"
    shell:
        """
        sourmash scripts fastmultigather --scaled 10_000 -k 51  {input.zip} {input.db} -o {output.gather}
        """

rule sourmash_tax_gtdb:
    input:
        gather=f"{outdir}/sourmash/{basename}-x-gtdb.fmg.csv",
        taxonomy_csv="/group/ctbrowngrp5/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv"
    output:
        tax=f"{outdir}/sourmash/{basename}-x-gtdb.classifications.csv"
    params:
        prefix=f"{outdir}/sourmash/{basename}-x-gtdb"
    shell:
        """
        sourmash tax genome --gather-csv {input.gather} \
                            --taxonomy-csv {input.taxonomy_csv} \
                            -o {params.prefix}
        """


rule sourmash_tax_ncbi:
    input:
        gather=f"{outdir}/sourmash/{basename}-x-ncbi.fmg.csv",
        taxonomy_csv="/group/ctbrowngrp5/sourmash-db/entire-2025-01-21/entire-2025-01-21.lineages.csv"
    output:
        tax=f"{outdir}/sourmash/{basename}-x-ncbi.classifications.csv"
    params:
        prefix=f"{outdir}/sourmash/{basename}-x-ncbi"
    shell:
        """
        sourmash tax genome --gather-csv {input.gather} \
                            --taxonomy-csv {input.taxonomy_csv} \
                            -o {params.prefix}
        """
# rule sendsketch:
#     input:
#         fasta="multi/bins/MEGAHIT-MetaBAT2-{sample}.{bin}.fa.gz"
#     output:
#         txt="results/sendsketch/{sample}.txt"
#     conda:
#         "envs/sendsketch.yml"
#     shell:
#         """
#         bbmap/sendsketch.sh in={input.fasta} refseq > {output.txt}
#         """

# rule gtdbtk_classify:
#     input:
#         fasta="bins/{sample}.fasta"
#     output:
#         summary="results/gtdbtk/{sample}/gtdbtk.bac120.summary.tsv"
#     params:
#         outdir=lambda wildcards: f"results/gtdbtk/{wildcards.sample}"
#     conda:
#         "envs/gtdbtk.yml"
#     shell:
#         """
#         gtdbtk classify_wf --genome_dir bins \
#                            --genome {input.fasta} \
#                            --out_dir {params.outdir} \
#                            --extension fasta \
#                            --cpus 4
#         """

# rule summarize_classification:
#     input:
#         sendsketch=expand("results/sendsketch/{sample}.txt", sample=BINS),
#         gtdbtk=expand("results/gtdbtk/{sample}/gtdbtk.bac120.summary.tsv", sample=BINS),
#         sourmash=expand("results/sourmash/{sample}.taxonomy.csv", sample=BINS)
#     output:
#         summary="results/summary/combined_taxonomy_summary.csv"
#     params:
#         bins=BINS
#     shell:
#         """
#         scripts/summarize_taxonomy.py \
#             --bins {params.bins} \
#             --sendsketch_dir results/sendsketch \
#             --gtdbtk_dir results/gtdbtk \
#             --sourmash_dir results/sourmash \
#             -o {output.summary}
#         """
