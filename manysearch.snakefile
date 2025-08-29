
outdir = "output.manysearch"
metagenome_mf = "br-magtest.wort.mf.csv",

search_genomes = {"GCF_000888735.1": "Acanthamoeba polyphaga mimivirus",
                  "GCF_001890705.1": "Aspergillus sydowii",
                  "GCF_003013715.1": "Candida auris",
                  "GCF_000165345.1": "Cryptosporidium parvum ",
                  "GCF_002217175.1": "Folsomia candida"}
# metagenomes = [
#     "ERR10162411", "ERR2819892", "SRR11734772", "SRR11734785", "SRR12217734", "SRR12959777", "SRR14842381", "SRR17231405", 
#     "SRR18691112", "SRR18691152", "SRR18691161", "SRR20285055", "SRR20950787", "SRR2243572", "SRR3458563", "SRR5098319", 
#     "SRR5190220", "SRR5190256", "SRR6144753", "SRR9016984", "ERR2185279", "ERR3333611", "SRR11734780", "SRR11734791", 
#     "SRR12217748", "SRR14141927", "SRR16797981", "SRR17238487", "SRR18691151", "SRR18691155", "SRR20285028", "SRR20950657", 
#     "SRR2105903", "SRR3458562", "SRR5098299", "SRR5190219", "SRR5190255", "SRR5831603", "SRR6144754", "SRR9016985",
# ]


rule all:
    input:
        f"{outdir}/search-genomes.sig.zip",
        f"{outdir}/search-genomes-x-brmetagenomes.manysearch.csv",
        # f"{outdir}/search-genomes-x-brmetagenomes.manysearch.checked.csv",
        f"{outdir}/bins-x-search-genomes.multisearch.csv",
        f"{outdir}/bins-x-search-genomes.multisearch.sc1000.csv",
        f"{outdir}/bins-x-search-genomes.multisearch.sc10.csv",
        f"{outdir}/bins-x-ncbi-entire.manysearch.csv",
        # f"{outdir}/bins-x-ncbi-euks.manysearch.sc1000.csv",
        #f"{outdir}/bins-x-ncbi-euks.prefetch.sc1000.csv",
        f"{outdir}/bins-x-ncbi-euks.multisearch.sc1000.csv",



rule write_directsketch_csv:
    output:
        csv=f"{outdir}/search-genomes.gbsketch.csv"
    run:
        with open(output.csv, 'w') as f:
            f.write("accession,name\n")
            for acc, species_name in search_genomes.items():
                name = f"{acc} {species_name}"
                f.write(f"{acc},{name}\n")


rule sourmash_sketch_search_genomes:
    input:
        csv=f"{outdir}/search-genomes.gbsketch.csv",
    output:
        sig=f"{outdir}/search-genomes.sig.zip",
    threads: 4
    benchmark: f"{outdir}/logs/sketch-search-genomes.benchmark"
    shell:
        """
        sourmash scripts gbsketch -p dna,k=21,k=31,k=51,scaled=1000 {input.csv} -o {output.sig}
        """

rule sourmash_manysearch_metagenomes:
    """
    search genomes against a set of metagenomes using sourmash manysearch
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sig.zip",
        metagenomes_mf = "br-magtest.wort.mf.csv",
    output:
        csv=f"{outdir}/search-genomes-x-brmetagenomes.manysearch.csv"
    threads: 4
    log: f"{outdir}/logs/search-genomes-x-brmetagenomes.manysearch.log"
    benchmark: f"{outdir}/logs/search-genomes-x-brmetagenomes.manysearch.benchmark"
    shell:
        """
        sourmash scripts manysearch {input.genomes_zip} {input.metagenomes_mf} -k 31 --scaled 1000 --output {output.csv} > {log} 2>&1
        """

rule write_manysketch_bins_csv:
    input:
        bins="multi_bins.txt",
    output:
        csv=f"{outdir}/bins.manysketch.csv"
    run:
        with open(output.csv, 'w') as f:
            f.write("name,genome_filename,protein_filename\n")
            # iterate over input fasta files and write to csv
            bin_list = [i.strip() for i in open(input.bins).readlines()]
            for bin_file in bin_list:
                sample_bin = os.path.basename(bin_file).replace("MEGAHIT-MetaBAT2-", "").replace(".fa", "")
                f.write(f"{sample_bin},{bin_file},\n")

rule manysketch_bins:
    input:
        csv=f"{outdir}/bins.manysketch.csv"
    output:
        sig=f"{outdir}/bins.sig.zip"
    threads: 4
    benchmark: f"{outdir}/logs/sketch-bins.benchmark"
    shell:
        """
        sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=100 {input.csv} -o {output.sig}
        """

rule manysketch_bins_sc10:
    input:
        csv=f"{outdir}/bins.manysketch.csv"
    output:
        sig=f"{outdir}/bins.sc10.sig.zip"
    threads: 4
    benchmark: f"{outdir}/logs/sketch-bins.benchmark"
    shell:
        """
        sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=10 {input.csv} -o {output.sig}
        """

# better - just do the orig sketching as higher res.
rule sourmash_sketch_search_genomes_highres:
    input:
        csv=f"{outdir}/search-genomes.gbsketch.csv",
    output:
        sig=f"{outdir}/search-genomes.sc10.sig.zip",
    threads: 4
    benchmark: f"{outdir}/logs/sketch-search-genomes.sc10.benchmark"
    shell:
        """
        sourmash scripts gbsketch -p dna,k=21,k=31,k=51,scaled=10 {input.csv} -o {output.sig}
        """

# now search bins against the search genomes
rule multisearch_bins_sc100:
    """
    search bins against the search genomes using sourmash manysearch
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sc10.sig.zip",
        bins_sig=f"{outdir}/bins.sig.zip",
    output:
        csv=f"{outdir}/bins-x-search-genomes.multisearch.csv"
    threads: 4
    log: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc100.log"
    benchmark: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc100.benchmark"
    shell:
        """
        sourmash scripts multisearch {input.bins_sig} {input.genomes_zip} -k 31 --scaled 100 --output {output.csv} --ani -m DNA --threshold 0.001 > {log} 2>&1
        """

rule multisearch_bins_sc10:
    """
    search bins against the search genomes using sourmash manysearch
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sc10.sig.zip",
        bins_sig=f"{outdir}/bins.sc10.sig.zip",
    output:
        csv=f"{outdir}/bins-x-search-genomes.multisearch.sc10.csv"
    threads: 4
    log: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc10.log"
    benchmark: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc10.benchmark"
    shell:
        """
        sourmash scripts multisearch {input.bins_sig} {input.genomes_zip} -k 31 --scaled 10 --output {output.csv} --ani -m DNA --threshold 0.01 > {log} 2>&1
        """

# do we really need the highres? try at sc1000
rule multisearch_bins_sc1000:
    """
    search bins against the search genomes using sourmash manysearch
    """
    input:
        genomes_zip=f"{outdir}/search-genomes.sig.zip",
        bins_sig=f"{outdir}/bins.sig.zip",
    output:
        csv=f"{outdir}/bins-x-search-genomes.multisearch.sc1000.csv"
    threads: 4
    log: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc1000.log"
    benchmark: f"{outdir}/logs/bins-x-search-genomes.multisearch.sc1000.benchmark"
    shell:
        """
        sourmash scripts multisearch {input.bins_sig} {input.genomes_zip} -k 31 --scaled 1000 --output {output.csv} --ani -m DNA > {log} 2>&1
        """


rule sourmash_manysearch_ncbi_entire:
    input:
        zip=f"{outdir}/bins.sig.zip",
        db= "/group/ctbrowngrp5/sourmash-db/entire-2025-01-21/entire-2025-01-21.k51.rocksdb"
    output:
        manysearch=f"{outdir}/bins-x-ncbi-entire.manysearch.csv"
    log: f"{outdir}/logs/bins-x-ncbi-entire.manysearch.log"
    benchmark: f"{outdir}/logs/bins-x-ncbi-entire.manysearch.benchmark"
    shell:
        """
        sourmash scripts manysearch --scaled 10_000 -k 51 {input.zip} {input.db} -o {output.manysearch} > {log} 2>&1
        """

rule sourmash_multisearch_ncbi_euk_zips:
    input:
        zip=f"{outdir}/bins.sig.zip",
        fungi="/group/ctbrowngrp5/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k51.zip",
        euk_other="/group/ctbrowngrp5/2025-genbank-eukaryotes/eukaryotes-other.sig.zip",
        bilateria_minus_verts="/group/ctbrowngrp5/2025-genbank-eukaryotes/bilateria-minus-vertebrates.sig.zip",
        db_txt="euk-dbs.txt"
    output:
        multisearch=f"{outdir}/bins-x-ncbi-euks.multisearch.sc1000.csv"
    log: f"{outdir}/logs/bins-x-ncbi-entire.multisearch.sc1000.log"
    benchmark: f"{outdir}/logs/bins-x-ncbi-euk.multisearch.sc1000.benchmark"
    shell:
        """
        sourmash scripts multisearch --scaled 1_000 -m DNA -k 51 --ani {input.zip} {input.db_txt} -o {output.multisearch} > {log} 2>&1
        """

rule sourmash_manysearch_ncbi_euk_zips:
    input:
        zip=f"{outdir}/bins.sig.zip",
#        fungi="/group/ctbrowngrp5/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k51.zip",
#        euk_other="/group/ctbrowngrp5/2025-genbank-eukaryotes/eukaryotes-other.sig.zip",
#        bilateria_minus_verts="/group/ctbrowngrp5/2025-genbank-eukaryotes/bilateria-minus-vertebrates.sig.zip",
        db_txt="euk-dbs.txt"
    output:
        manysearch=f"{outdir}/bins-x-ncbi-euks.manysearch.sc1000.csv"
    log: f"{outdir}/logs/bins-x-ncbi-entire.manysearch.sc1000.log"
    benchmark: f"{outdir}/logs/bins-x-ncbi-euk.manysearch.sc1000.benchmark"
    shell:
        """
        sourmash scripts manysearch --scaled 1_000 -k 51 {input.zip} {input.db_txt} -o {output.manysearch} > {log} 2>&1
        """
        #sourmash scripts manysearch --scaled 1_000 -k 51 {input.zip} {input.fungi} {input.euk_other} {input.bilateria_minus_verts} -o {output.manysearch} > {log} 2>&1

# now aggregate and assess the results