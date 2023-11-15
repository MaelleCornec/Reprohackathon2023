# Chargement du fichier de configuration config.yaml
configfile:
    "config.yaml"

# Snakemake définit par défaut la première règle comme la cible.
# A l'état actuel, je n'ai pas encore inclus les script python dans le snakefile,
# Donc ma cible est le fichier count.
# Comme ça, quand on exécute la commande snakemake --cores, 
# Il exécutera récursivement toutes les règles nécessaires pour obtenir le fichier count
rule all:
    input:
        "data_comptage/counts.txt"

# Téléchargement du génome de référence
# dans un dossier reference
####### Ce serait une bonne idée de mettre le lien dans un autre fichier 
####### Pour pouvoir plus aptement utiliser ce scrip ailleurs
rule genome:
    output:
        "genome/reference.fasta"
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """

# Création de l'index sur le génome de référence
rule indexing:
    input:
        "genome/reference.fasta"
    output:
        "index/indexation"
    container:
        "docker://suzannegtx/bowtie:0.12.7"
    shell: 
        "bowtie-build {input} {output}"

def input_fastq(wildcards):
    return config["samples"][wildcards.sample]

# Découpe des morceaux d'échantillons qui sont trop incertains
# et élimination des échantillons de moins de 25 nucléotides
rule trimming:
    input:
        "data40000/{sample}.fastq"
    output:
        "data_trim/{sample}_trimmed.fq"
    container:
        "docker://suzannegtx/trim-galore:0.6.4"
    shell:
        """
        trim_galore -q 20 --phred33 --length 25 {input} -O data_trim
        """

# Cartographie des échantillons grâce à l'index fait sur le génome de référence
rule mapping:
    input:
        ind="index/indexation",
        data="data_trim/{sample}_trimmed.fq"      
    output:
        "data_map/{sample}.bam"
    container:
        "docker://suzannegtx/trim-galore:0.6.4"
    shell:
        """
        bowtie -p 4 -S -x {input.ind} {input.data} | samtools sort -@ 4 -o {output}
        samtools index {output}
        """

# Téléchargement des annotations du génome de référence
# dans un dossier reference
####### Ce serait une bonne idée de mettre le lien dans un autre fichier 
####### Pour pouvoir plus aptement utiliser ce scrip ailleurs
rule annotation_gen:
    output:
        "genome/reference.gff"
    shell:
        """
        wget -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
        """

# Dénombrement des motifs présents sur les échantillons
# grâce aux annotations du génome de référence
rule counting:
    input:
        ech=expand("data_map/{sample}.bam", echantillon=config["samples"]),
        ref="genome/reference.gff"
    output:
        "data_comptage/counts.txt"
    container:
        "docker://suzannegtx/subreads-featurecounts:1.4.6-p3"
    shell:
        "featureCounts --extraAttributes Name -t gene -g -s 1 ID -F GTF -T 4 -a {input.ref} -o {output} {input.ech}"

# Analyse des échantillons
# en connaissant le nom des gènes
#rule analysis:
#    input:
#        "data40000compt/counts.txt"
#    output:
#        "je ne sais vraiment plus"
#    script:
#        "analysis.py"