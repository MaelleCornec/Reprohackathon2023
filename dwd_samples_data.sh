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
SRA_ACCESSION_1="SRR10379721"
SRA_ACCESSION_2="SRR10379722"
SRA_ACCESSION_3="SRR10379723"
SRA_ACCESSION_4="SRR10379724"
SRA_ACCESSION_5="SRR10379725"
SRA_ACCESSION_6="SRR10379726"

# Répertoire de sortie
OUTPUT_DIRECTORY="/home/ubuntu/data_echantillons"

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

# Renommer les fichiers FASTQ
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_1".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_1}.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_2".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_2}.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_3".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_3}.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_4".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_4}.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_5".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_5}.fastq"
mv "$OUTPUT_DIRECTORY/$SRA_ACCESSION_6".fastq "$OUTPUT_DIRECTORY/${SRA_ACCESSION_6}.fastq"