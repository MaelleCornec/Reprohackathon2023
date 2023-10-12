#!/bin/bash


# Vérification de la version installée
fastq-dump --version

cd data_echantillons

# Numéro d'accès SRA
SRA_ACCESSION_6="SRR10379726"  #control replicate 3

# Répertoire de sortie
OUTPUT_DIRECTORY="/home/ubuntu/Reprohackathon2023/data_echantillons"

# Nom de fichier de sortie personnalisé (sans extension .fastq)

NOM_FICHIER_6="control_3"


# Téléchargement des fichiers SRA avec prefetch

prefetch "$SRA_ACCESSION_6"

# Conversion des fichiers SRA en FASTQ avec renommage (sans compression)
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_6" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_6"

# Renommer les fichiers FASTQ (paired-end)
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_1"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_1}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_1"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_1}_R2.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_2"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_2}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_2"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_2}_R2.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_3"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_3}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_3"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_3}_R2.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_4"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_4}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_4"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_4}_R2.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_5"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_5}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_5"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_5}_R2.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_6"_1.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_6}_R1.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_6"_2.fastq "$OUTPUT_DIRECTORY/${NOM_FICHIER_6}_R2.fastq"