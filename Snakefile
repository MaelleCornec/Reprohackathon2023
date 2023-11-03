# Référencement du nom des 6 fichiers à étudier
ECHANTILLONS = ['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726']

# Définition de la wildcard echantillon
rule all:
    input:
        expand("data40000/{echantillon}.fastq", echantillon=ECHANTILLONS)

# Téléchargement du génome de référence
# dans un dossier reference
rule genome:
    output:
        "reference/genome.fasta"
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """

# Création de l'index sur le génome de référence
rule indexing:
    input:
        "reference/genome.fasta"
    output:
        "index/indexation"
    shell: 
        "bowtie-build {input} {output}"

# Découpe des morceaux d'échantillons qui sont trop incertains
# et élimination des échantillons de moins de 25 nucléotides
rule trimming:
    input:
        expand("data40000/{echantillon}.fastq", echantillon=ECHANTILLONS)
    output:
        expand("{echantillon}_trimmed.fq {echantillon}.fastq_trimming_report.txt", echantillon=ECHANTILLONS)
    shell:
        """
        mkdir data40000trim
        trim_galore -q 20 --phred33 --length 25 {input}
        mv {output} data40000trim
        """

# Cartographie des échantillons grâce à l'index fait sur le génome de référence
rule mapping:
    input:
        expand("data40000trim/{echantillon}_trimmed.fq", echantillon=ECHANTILLONS)      
    output:
        expand("data40000map/{echantillon}.bam", echantillon=ECHANTILLONS)
    shell:
        """
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379721_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379721.bam
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379722_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379722.bam
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379723_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379723.bam
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379724_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379724.bam
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379725_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379725.bam
        bowtie -p 4 -S -x index/indexation data40000trim/SRR10379726_trimmed.fq | samtools sort -@ 4 -o data40000map/SRR10379726.bam
        samtools index data40000map/SRR10379721.bam
        samtools index data40000map/SRR10379722.bam
        samtools index data40000map/SRR10379723.bam
        samtools index data40000map/SRR10379724.bam
        samtools index data40000map/SRR10379725.bam
        samtools index data40000map/SRR10379726.bam
        """

# Téléchargement des annotations du génome de référence
# dans un dossier reference
rule annotation_gen:
    output:
        "reference/genome.gff"
    shell:
        """
        wget -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
        """

# Dénombrement des motifs présents sur les échantillons
# grâce aux annotations du génome de référence
rule counting:
    input:
        ech=expand("data40000map/{echantillon}.bam", echantillon=ECHANTILLONS),
        ref="reference/genome.gff"
    output:
        "data40000compt/counts.txt"
    shell:
        "featureCounts --extraAttributes Name -t gene -g ID -F GTF -T 4 -a {input.ref} -o {output} {input.ech}"

# Analyse des échantillons
# en connaissant le nom des gènes
#rule analysis:
#    input:
#        "data40000compt/counts.txt"
#    output:
#        "je ne sais vraiment plus"
#    script:
#        "analysis.py"