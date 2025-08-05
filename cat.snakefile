# Snakemake workflow for BAT (BIN Annotation Tool) annotation against NR database
# This workflow processes multiple metagenomes, annotating bins using CAT_pack and summarizing results


metagenomes = [
    "ERR10162411", "ERR2819892", "SRR11734772", "SRR11734785", "SRR12217734", "SRR12959777", "SRR14842381", "SRR17231405", 
    "SRR18691112", "SRR18691152", "SRR18691161", "SRR20285055", "SRR20950787", "SRR2243572", "SRR3458563", "SRR5098319", 
    "SRR5190220", "SRR5190256", "SRR6144753", "SRR9016984", "ERR2185279", "ERR3333611", "SRR11734780", "SRR11734791", 
    "SRR12217748", "SRR14141927", "SRR16797981", "SRR17238487", "SRR18691151", "SRR18691155", "SRR20285028", "SRR20950657", 
    "SRR2105903", "SRR3458562", "SRR5098299", "SRR5190219", "SRR5190255", "SRR5831603", "SRR6144754", "SRR9016985",
]
# metagenomes = ["SRR11734772"]
# metagenomes = ["ERR10162411"]

outdir = "output.BAT"
thisdir = "/home/ntpierce/2025-br-mag"

rule all:
    input: 
        expand(f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.summary.txt", metagenome=metagenomes),
        expand(f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.taxnames.txt", metagenome=metagenomes)


rule bat_bin_annotation:
    input:
        bins=f"{thisdir}/multi/bins/{{metagenome}}"
    output:
        # Specify the output file(s) here if applicable
        fasta=f"{outdir}/{{metagenome}}/out.BAT.concatenated.fasta",
        log=f"{outdir}/{{metagenome}}/out.BAT.log",
        prodigal=f"{outdir}/{{metagenome}}/out.BAT.concatenated.predicted_proteins.faa",
        prodigal_gff=f"{outdir}/{{metagenome}}/out.BAT.concatenated.predicted_proteins.gff",
        diamond=f"{outdir}/{{metagenome}}/out.BAT.concatenated.alignment.diamond",
        orf2lca=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.txt",
        classif=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.txt",
    params:
        batdir=f"{outdir}/{{metagenome}}",
        db=f"{thisdir}/20241212_CAT_nr_website/db/",
        tax=f"{thisdir}/20241212_CAT_nr_website/tax/",
        suffix=".fa",
        # BAT output files
        # fasta="out.BAT.concatenated.fasta",
        # log="out.BAT.log",
        # prodigal="out.BAT.concatenated.predicted_proteins.faa",
        # prodigal_gff="out.BAT.concatenated.predicted_proteins.gff",
        # diamond="out.BAT.concatenated.alignment.diamond",
        # orf2lca="out.BAT.ORF2LCA.txt",
        # classif="out.BAT.bin2classification.txt",
    benchmark: f"{outdir}/benchmarks/{{metagenome}}.BAT.benchmark"
    threads: 20
    resources:
        mem_mb=20000,
        time="3-00:00:00",
        slurm_partition="med2",
    shell:
        """
        cd {params.batdir}
        CAT_pack bins -b {input.bins} -d {params.db} -t {params.tax} -s {params.suffix:q}
        """
        #   mv {params.fasta} {output.fasta}
        # mv {params.log} {output.log}
        # mv {params.prodigal} {output.prodigal}
        # mv {params.prodigal_gff} {output.prodigal_gff}
        # mv {params.diamond} {output.diamond}
        # mv {params.orf2lca} {output.orf2lca}
        # mv {params.classif} {output.classif}


rule add_taxnames:
    input:
        tax="20241212_CAT_nr_website/tax/",
        orf2lca=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.txt",
        classif=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.txt",
    output:
        taxnames=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.taxnames.txt",
        classif_taxnames=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.taxnames.txt",
    benchmark: f"{outdir}/benchmarks/{{metagenome}}.add-taxnames.benchmark"
    shell:
        """
        CAT_pack add_names -i {input.orf2lca} -t {input.tax} -o {output.taxnames}
        CAT_pack add_names -i {input.classif} -t {input.tax} -o {output.classif_taxnames} --only_official
        """

rule summarize_classification:
    input:
        classif_taxnames=f"{outdir}/{{metagenome}}/out.BAT.bin2classification.taxnames.txt",
    output:
        summary=f"{outdir}/{{metagenome}}/out.BAT.ORF2LCA.summary.txt",
    benchmark: f"{outdir}/benchmarks/{{metagenome}}.summarize-classification.benchmark"
    shell:
        """
        CAT_pack summarise -i {input.classif_taxnames} -o {output.summary}
        """
