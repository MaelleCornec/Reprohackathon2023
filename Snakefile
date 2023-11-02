rule indexing:
    input:
        "reference.fasta"
    output:
        "indexation"
    shell: 
        "bowtie-build genome.fastq indexation"