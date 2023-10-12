# ce qu'il faut sortir à la fin
# je ne sais pas
# on n'est pas encore arrivé jusque-là

echantillons = ['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726']

rule all:
    input:
        expand("....", echantillon=echantillons)

# téléchargement du génome de référence
# on peut utiliser le même conteneur que pour le téléchargement des échantillons, non ?

rule downloadingenome:
    output:
        "reference/genome.fastq"
    container:
        "fastqc:v1"
    threads:1
    script:
        "downloadingen.sh"

# dans downloadingen.sh :
#$ wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

# téléchargement des annotations du génome de référence
# je ne sais pas du tout où trouver ça ?
# et quelle est la forme de ce fichier d'annotation ?

rule downloadingannot:
    output:
        "annotation.txt"
    script:
        "downloadingannot.sh"

# dans downloadingannot.sh :
# $ wget -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&→report=gff3&id=CP000253.1"

# téléchargement des noms de gènes du génome de référence et des chemins
# je ne sais pas du tout où trouver ça ?
# et quelle est la forme de ce fichier d ?

rule downloadingnoms:
    output:
        "noms.txt"
    script:
        "downloadingnoms.py"

# création de l'index sur le génome de référence

rule indexing:
    input:
        "genome.fastq"
    output:
        "indexGenome ???"
    container:
        "bowtie:v1"
    script:
        "indexing.sh"

# dans indexing.sh :
# $ bowtie-build genome.fastq index

# téléchargement des données 6 fichiers fastq
# je veux que mes 6 échantillons soient dans un dossier nommé 'echantillons'
# mais est-ce que c'est vraiment pratique de retélécharger ces fichiers à chaque fois ? Il sont très volumineux...
# ah, mais snakemake ne devrait pas me les retélécharger s'ils sont déjà là, donc c'est pas grave
# je continue

rule downloadingechantillons:
    input:
        "{echantillon}""
    output: 
        "echantillons/{echantillon}.fastq"
    container:
        "fastqc:v1"
    threads:4
    script: 
        "downloadingech.sh"

# dans downloadingech.sh (à développer, voir script de Suzanne):
# $ fasterq-dump --progress {echantillon}
# $ gzip *.fastq

# Trimming : coupe des morceaux d'échantillons qui sont trop incertains
# et élimination des échantillons de moins de 25 nucléotides

rule trimming:
    input:
        "liste des echantillons : echantillonInit/{echantillon}.fastq"
    output:
        "liste des échantillons moins ceux qui sont trop courts/peu certains : echantillonCoupe/{echantillon}.fastq"
    container:
        "trimGalore:v1"
    script:
        "trimming.sh"

# dans trimming.sh (à développer):
# $ trim_galore -q 20 --phred33 --length 25 echantillonInit/{echantillon}.fastq

# Cartographie des échantillons grâce à l'index fait sur le génôme de référence

rule mapping:
    input:
        "échantillons coupés précédemment : echantillonCoupe/{echantillon}.fastq"
        "indexGenome"
    output:
        "échantillons cartographiés : echantillonCarto/{echantillon}.fastq"
    container:
        "bowtie:v1"
    script:
        "mapping.sh"

# dans mapping.sh (à développer et à comprendre):
# $ bowtie -p <#cpus> -S indexGenome <(gunzip -c <GZIPED FASTQ FILE>) | \samtools sort -@ <#CPUS> > <NAME>.bam
# $ samtools index <NAME>.bam

# dénombrage des motifs présents sur les échantillons
# grâce aux annotations du génome de référence

rule counting:
    input:
        "echantillonCarto/{echantillon}.fastq"
        "annotation.txt"
    output:
        "echantillonAnnot/{echantillon}.fastq ???"
    container:
        "featureCount:v1"
    script:
        "counting.sh"

# dans counting.sh
# $ featureCounts --extraAttributes Name -t gene -g ID -F GTF -T <#CPUS> -a <GFF> -o counts.txt <BAM FILES>

# analyse des échantillons
# en connaissant le nom des gènes

rule analysis:
    input:
        "echantillonAnnot/{echantillon}.fastq ???"
        "nom.txt"
    output:
        "je ne sais vraiment plus"
    container:
        "deSeq2:v1"
    script:
        "analysis.py"

# en vrai, je pourrais probablement mettre tous les scripts directement dans les étapes
# en utilisant : 
# shell: " ... "
