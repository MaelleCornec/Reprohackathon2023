# des téléchargements ont déjà été faits :
# - les 6 échantillons d'intérêt (6 fichiers fastq)
# - le génome de référence (reference.fasta)
# - les annotations du génome de  référence (reference.gff)
# - les noms des gènes du génome de référence (on ne les a pas encore !)

# des containers ont été créés :
# - trim_galore
# - bowtie
# - feature_count

echantillons = ['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726']

rule all:
    input:
        "reference.fasta"
        "reference.gff"
        expand("echantillonstests/{echantillon}.fastq", echantillon=echantillons)

# création de l'index sur le génome de référence

rule indexing:
    input:
        "reference.fasta"
    output:
        "indexation.txt"
    container:
        "bowtie"
    shell: 
        "bowtie-build genome.fastq indexation.txt"

# Trimming : coupe des morceaux d'échantillons qui sont trop incertains
# et élimination des échantillons de moins de 25 nucléotides

rule trimming:
    input:
        "{echantillon}.fastq"
    output:
        "echantillonstests/{echantillon}.fastq"
    container:
        "trim_galore"
    shell:
        "trim_galore -q 20 --phred33 --length 25 echantillonstests/{echantillon}.fastq"

# Cartographie des échantillons grâce à l'index fait sur le génôme de référence

rule mapping:
    input:
        "{echantillon}.fastq"
        "indexation.txt"
    output:
        "carto/{echantillon}.bam"
    container:
        "bowtie"
    shell:
        """
        bowtie -p 4 -S indexation.txt {echantillon}.fastq | \samtools sort -@ 4 carto/{echantillon}.bam
        samtools index carto/{echantillon}.bam
        """

# dénombrage des motifs présents sur les échantillons
# grâce aux annotations du génome de référence

rule counting:
    input:
        "carto/{echantillon}.bam"
        "reference.gff"
    output:
        "compt/counts.txt"
    container:
        "feature_count"
    shell:
        "featureCounts --extraAttributes Name -t gene -g ID -F GTF -T 4 -a reference.gff -o counts.txt {echantillon}.bam

# analyse des échantillons
# en connaissant le nom des gènes

rule analysis:
    input:
        "compt/counts.txt"
    output:
        "je ne sais vraiment plus"
    container:
        "de_seq2"
    script:
        "analysis.py"

##########################
# Dénombrement des motifs présents sur les échantillons
# grâce aux annotations du génome de référence
rule counting:
    input:
        "data40000map/{echantillon}.bam",
        "reference.gff"
    output:
        "data40000compt/counts.txt"
    shell:
        "featureCounts --extraAttributes Name -t gene -g ID -F GTF -T 4 -a reference.gff -o counts.txt {echantillon}.bam

# analyse des échantillons
# en connaissant le nom des gènes

rule analysis:
    input:
        "compt/counts.txt"
    output:
        "je ne sais vraiment plus"
    script:
        "analysis.py"