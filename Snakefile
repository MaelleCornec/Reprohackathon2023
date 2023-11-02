# Référencement du nom des 6 fichiers à étudier
ECHANTILLONS = ['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726']

# Définition de la wildcard echantillon
rule all:
    input:
        expand("data40000/{echantillon}.fastq", echantillon=ECHANTILLONS)

# Création de l'index sur le génome de référence
rule indexing:
    input:
        "reference.fasta"
    output:
        "index/indexation"
    shell: 
        "bowtie-build {input} {output}"

# Découpe des morceaux d'échantillons qui sont trop incertains
# et élimination des échantillons de moins de 25 nucléotides
rule trimming:
    input:
        "data40000/{echantillon}.fastq"
    output:
        "data40000trim/{echantillon}_trimmed.fq"
    shell:
        "trim_galore -q 20 --phred33 --length 25 {input} {output}"

# Cartographie des échantillons grâce à l'index fait sur le génome de référence
rule mapping:
    input:
        "index/*",
        "data40000trim/{echantillon}_trimmed.fq"
    output:
        "data40000map/{echantillon}.bam"
    shell:
        """
        bowtie -p 4 -S {input} | \samtools sort -@ 4 {output}
        samtools index {output}
        """

