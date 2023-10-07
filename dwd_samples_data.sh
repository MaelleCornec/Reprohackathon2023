#!/bin/bash

#Mise à jour des informations sur les paquets
sudo apt-get update

# Installation de SRA Toolkit
sudo apt-get install sra-toolkit

# Vérification de la version installée
fastq-dump --version

mkdir data_echantillons
cd data_echantillons

# Numéro d'accès SRA
SRA_ACCESSION_1="SRR10379721"  #persister replicate 1
SRA_ACCESSION_2="SRR10379722"  #perister replicate 2
SRA_ACCESSION_3="SRR10379723"  #persister replicate 3
SRA_ACCESSION_4="SRR10379724"  #control replicate 1
SRA_ACCESSION_5="SRR10379725"  #control replicate 2
SRA_ACCESSION_6="SRR10379726"  #control replicate 3

# Répertoire de sortie
OUTPUT_DIRECTORY="/home/ubuntu/Reprohackathon2023/data_echantillons"

# Nom de fichier de sortie personnalisé (sans extension .fastq)
NOM_FICHIER_1="persister_replicate_1"
NOM_FICHIER_2="persister_replicate_2"
NOM_FICHIER_3="persister_replicate_3"
NOM_FICHIER_4="control_1"
NOM_FICHIER_5="control_2"
NOM_FICHIER_6="control_3"


# Téléchargement des fichiers SRA avec prefetch
prefetch "$SRA_ACCESSION_1"
prefetch "$SRA_ACCESSION_2"
prefetch "$SRA_ACCESSION_3"
prefetch "$SRA_ACCESSION_4"
prefetch "$SRA_ACCESSION_5"
prefetch "$SRA_ACCESSION_6"

# Conversion des fichiers SRA en FASTQ avec renommage (sans compression)
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_1" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_1"
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_2" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_2"
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_3" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_3"
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_4" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_4"
fastq-dump --outdir "$OUTPUT_DIRECTORY" --split-files --origfmt --defline-seq '@$sn/$ri' --defline-qual '+' "$SRA_ACCESSION_5" -O "$OUTPUT_DIRECTORY/$NOM_FICHIER_5"
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