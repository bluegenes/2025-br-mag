# 2025-br-mag


Check Metagenome Samples and Metagenome Assembled Genome sequences used as part of branchwater-web ms.

1. Install and activate conda environment

```
mamba env create -f environment.yml
mamba activate 2025-br-mag
```

2. Run sourmash analyses

Download sourmash NCBI databases, if desired:

*links to be added*

```
snakemake -s sourmash.snakefile -n
```
> remove `-n` to run
> This analysis will run faster with additional cores, e.g. 4-30.

3. Run BAT annotation for bins

Download bin files:

*links to be added*

Run BAT analysis:

```
snakemake -s BAT.snakefile -n
```