# convert a bunch of .sig files into .sig.zip files and also produce .mf.csv files.
import os

metagenomes = [
    "ERR10162411", "ERR2819892", "SRR11734772", "SRR11734785", "SRR12217734", "SRR12959777", "SRR14842381", "SRR17231405", 
    "SRR18691112", "SRR18691152", "SRR18691161", "SRR20285055", "SRR20950787", "SRR2243572", "SRR3458563", "SRR5098319", 
    "SRR5190220", "SRR5190256", "SRR6144753", "SRR9016984", "ERR2185279", "ERR3333611", "SRR11734780", "SRR11734791", 
    "SRR12217748", "SRR14141927", "SRR16797981", "SRR17238487", "SRR18691151", "SRR18691155", "SRR20285028", "SRR20950657", 
    "SRR2105903", "SRR3458562", "SRR5098299", "SRR5190219", "SRR5190255", "SRR5831603", "SRR6144754", "SRR9016985",
    "ERR1992808", "ERR5083269", "SRR22535178", # these three did not produce no bins
]

sig_set_basename = "br-magtest"

rule all:
    input:
        expand("metagenomes/{acc}.sig.zip", acc=metagenomes),
        expand("metagenomes/{acc}.mf.csv", acc=metagenomes),
        f"{sig_set_basename}.wort.mf.csv",
    

rule make_sig_zip:
    input: "/group/ctbrowngrp/irber/data/wort-data/wort-sra/sigs/{acc}.sig"
    output: "metagenomes/{acc}.sig.zip"
    resources:
        slurm_partition="low2",
        mem_mb=lambda wildcards, attempt: attempt * 3000, # 3 GB per attempt
        runtime=60,
    threads: 1
    shell:
        """
        sourmash sig cat {input} -o {output}
        """

rule make_mf_csv:
    input: "metagenomes/{acc}.sig.zip"
    output: "metagenomes/{acc}.mf.csv"
    resources:
        slurm_partition="low2",
        mem_mb=lambda wildcards, attempt: attempt * 3000, # 3 GB per attempt
        runtime=240,
    shell:
        """
        sourmash sig collect {input} -o {output} -F csv --abspath
        """

rule aggregate_mf_csv:
    input: expand("metagenomes/{acc}.mf.csv", acc=metagenomes)
    output: f"{sig_set_basename}.wort.mf.csv"
    resources:
        slurm_partition="low2",
        mem_mb=lambda wildcards, attempt: attempt * 5000, # 5 GB per attempt
        runtime=240,
    shell:
        """
        sourmash sig collect {input} -o {output} -F csv --relpath
        """
