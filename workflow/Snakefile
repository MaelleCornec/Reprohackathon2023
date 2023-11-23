# Chargement du fichier de configuration config.yaml
configfile:
    "config.yaml"

# Snakemake définit par défaut la première règle comme la cible.
# A l'état actuel, je n'ai pas encore inclus les script python dans le snakefile, donc ma cible est le fichier count.
# Comme ça, quand on exécute la commande snakemake --cores, 
# Il exécutera récursivement toutes les règles nécessaires pour obtenir le fichier count
rule all:
    input:
        "data_comptage/counts.txt"

# Téléchargement des 6 échantillons à analyser
rule downloading:
    output:
        "data/{sample}.fastq"
    message:
        "Téléchargement de {wildcards.sample}"
    container:
        "docker://suzannegtx/sra-toolkit:2.10.0"
    log:
        "logs/downloading/{sample}.log"
    shell:
        """
        (fasterq-dump --threads 8 --progress {wildcards.sample} -O data) 2> {log}
        """

# Découpe des morceaux d'échantillons qui sont trop incertains et élimination des échantillons de moins de 25 nucléotides
rule trimming:
    input:
        "data/{sample}.fastq"
    output:
        "data_trim/{sample}_trimmed.fq"
    message:
        "Nettoyage de {wildcards.sample}"
    container:
        "docker://suzannegtx/trimgalore-cutadapt:0.6.4"
    log:
        "logs/trimming/{sample}.log"
    shell:
        """
        (trim_galore -q 20 --phred33 --length 25 {input} -O data_trim) 2> {log}
        """

# Téléchargement du génome de référence dans un dossier reference
rule genome:
    output:
        "genome/reference.fasta"
    params:
        config["linkfasta"]
    message:
        "Téléchargement du génome de référence"
    log:
        "logs/downloading/reference/fasta.log"
    shell:
        """
        (wget -q -O {output} "{params}") 2> {log}
        """

# Création de l'index sur le génome de référence
rule indexing:
    input:
        "genome/reference.fasta"
    output:
        "index/indexation"
    message:
        "Indexation du génome de référence"
    container:
        "docker://suzannegtx/bowtie:0.12.7"
    log:
        "logs/indexation.log"
    shell: 
        """
        (bowtie-build {input} {output}) 2> {log}
        """

# Création des fichiers sam
rule samming:
    input:
        "data_trim/{sample}_trimmed.fq"
    output:
        "data_sam/{sample}.sam"
    message:
        "Première étape du mapping de {wildcards.sample}"
    container:
        "docker://suzannegtx/bowtie:0.12.7"
    log:
        "logs/samming/{sample}.log"
    shell:
        """
        (bowtie -S -x index/indexation {input} > {output}) 2> {log}
        """

# Création des fichiers bam
rule bamming:
    input:
        "data_sam/{sample}.sam"
    output:
        "data_bam/{sample}.bam"
    message:
        "Deuwième étape du mapping de {wildcards.sample}"
    container:
        "docker://pegi3s/samtools_bcftools:1.9"
    log:
        "logs/bamming/{sample}.log"
    shell:
        """
        (samtools sort -O bam -o {output} {input}) 2> {log}
        """

# Cartographie des échantillons grâce à l'index fait sur le génome de référence
rule mapping:
    input:
        "data_bam/{sample}.bam"      
    output:
        "data_bam/{sample}.bam.bai"
    message:
        "Troisième étape du mapping de {wildcards.sample}"
    container:
        "docker://pegi3s/samtools_bcftools:1.9"
    log:
        "logs/mapping/{sample}.log"
    shell:
        """
        (samtools index {input}) 2> {log}
        """

# Téléchargement des annotations du génome de référence dans un dossier reference
rule annotation_gen:
    output:
        "genome/reference.gff"
    params:
        config["linkgff"]
    message:
        "Téléchargement des annotations du génome de référence"
    log:
        "logs/downloading/reference/gff.log" 
    shell:
        """
        (wget -O {output} "{params}") 2> {log}
        """

# Dénombrement des motifs présents sur les échantillons grâce aux annotations du génome de référence
rule counting:
    input:
        ech=expand("data_map/{sample}.bam", sample=config["samples"]),
        ref="genome/reference.gff"
    output:
        "data_comptage/counts.txt"
    message:
        "Dénombrement des motifs"
    container:
        "docker://suzannegtx/subreads-featurecounts:1.4.6-p3"
    log:
        "logs/counting.log"
    shell:
        """
        (featureCounts --extraAttributes Name -t gene -g ID -s 1 -F GTF -T 4 -a {input.ref} -o {output} {input.ech}) 2> {log}
        """

# Analyse des échantillons en connaissant le nom des gènes
#rule analysis:
#    input:
#        "data_comptage/counts.txt"
#    output:
#        "je ne sais vraiment plus"
#    script:
#        "analysis.py"
