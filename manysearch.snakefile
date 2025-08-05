
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

rule sourmash_manysearch:
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

