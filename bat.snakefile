# Snakemake workflow for BAT (BIN Annotation Tool) annotation against NR database
# This workflow processes multiple metagenomes, annotating bins using CAT_pack and summarizing results

## Instructions to run:
# 1. Make sure you're in the `2025-br-mag` directory and have the `2025-br-mag` conda environment activated.

# 2. Download a pre-prepared NR database (~290G). 

# You can download these from: https://tbb.bio.uu.nl/tina/CAT_pack_prepare/
# Here are the commands I used:

# wget https://tbb.bio.uu.nl/tina/CAT_pack_prepare/20241212_CAT_nr.tar.gz
# tar -xvzf 20241212_CAT_nr.tar.gz

# 3. Update the following database information paths:

database_dir = "20241212_CAT_nr_website/db/"
tax_dir = "20241212_CAT_nr_website/tax/"

# 4. Download the metagenome bins you want to annotate.

## TO DO: ADD DOWNLOAD COMMANDS FOR THE METAGENOME BINS ##

# 5. Modify `thisdir` to reflect your path to the `2025-br-mag` directory
# We're using this to ensure full paths to all files and databases, as we'll be running BAT from a different directory.
thisdir = "/home/ntpierce/2025-br-mag"

# 6. Execute the workflow with:
#    snakemake --use-conda --cores 20 --rerun-incomplete --keep-going --latency-wait 60 --printshellcmds

outdir = "output.BAT"

metagenomes = [
    "ERR10162411", "ERR2819892", "SRR11734772", "SRR11734785", "SRR12217734", "SRR12959777", "SRR14842381", "SRR17231405", 
    "SRR18691112", "SRR18691152", "SRR18691161", "SRR20285055", "SRR20950787", "SRR2243572", "SRR3458563", "SRR5098319", 
    "SRR5190220", "SRR5190256", "SRR6144753", "SRR9016984", "ERR2185279", "ERR3333611", "SRR11734780", "SRR11734791", 
    "SRR12217748", "SRR14141927", "SRR16797981", "SRR17238487", "SRR18691151", "SRR18691155", "SRR20285028", "SRR20950657", 
    "SRR2105903", "SRR3458562", "SRR5098299", "SRR5190219", "SRR5190255", "SRR5831603", "SRR6144754", "SRR9016985",
]

rule all:
    input: 
        expand(f"{outdir}/{{metagenome}}/out.BAT.bin2classification.taxnames.txt", metagenome=metagenomes),
        expand(f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.taxnames.txt", metagenome=metagenomes),
        # expand(f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.summary.txt", metagenome=metagenomes),

rule bat_bin_annotation:
    input:
        bins=f"{thisdir}/multi/bins/{{metagenome}}"
    output:
        fasta=f"{outdir}/{{metagenome}}/out.BAT.concatenated.fasta",
        log=f"{outdir}/{{metagenome}}/out.BAT.log",
        prodigal=f"{outdir}/{{metagenome}}/out.BAT.concatenated.predicted_proteins.faa",
        prodigal_gff=f"{outdir}/{{metagenome}}/out.BAT.concatenated.predicted_proteins.gff",
        diamond=f"{outdir}/{{metagenome}}/out.BAT.concatenated.alignment.diamond",
        orf2lca=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.txt",
        classif=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.txt",
    params:
        batdir=f"{outdir}/{{metagenome}}",
        db=f"{thisdir}/{database_dir}",
        tax=f"{thisdir}/{tax_dir}",
        suffix=".fa",
    benchmark: f"{outdir}/benchmarks/{{metagenome}}.BAT.benchmark"
    threads: 20
    conda: "bat-env.yml"
    resources:
        mem_mb=20000,
        time="3-00:00:00",
        slurm_partition="med2",
    shell:
        """
        cd {params.batdir}
        CAT_pack bins -b {input.bins} -d {params.db} -t {params.tax} -s {params.suffix:q}
        """

rule add_taxnames:
    input:
        tax=tax_dir,
        orf2lca=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.txt",
        classif=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.txt",
    output:
        orf_taxnames=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.taxnames.txt",
        classif_taxnames=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.taxnames.txt",
    benchmark: f"{outdir}/benchmarks/{{metagenome}}.add-taxnames.benchmark"
    conda: "bat-env.yml"
    shell:
        """
        CAT_pack add_names -i {input.classif} -t {input.tax} -o {output.classif_taxnames} --only_official
        CAT_pack add_names -i {input.orf2lca} -t {input.tax} -o {output.orf_taxnames}
        """

# rule summarize_classification:
#     input:
#         classif_taxnames=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.taxnames.txt",
#     output:
#         summary=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.summary.txt",
#     benchmark: f"{outdir}/benchmarks/{{metagenome}}.summarize-classification.benchmark"
#     conda: "bat-env.yml"
#     shell:
#         """
#         CAT_pack summarise -i {input.classif_taxnames} -o {output.summary}
#         """
