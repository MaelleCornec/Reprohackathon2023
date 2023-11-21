# Spécification de l'image de base Ubuntu
FROM ubuntu:latest

# Mise à jour des packages et installation des dépendances nécessaires
RUN apt update && \
    apt install -y wget unzip

# Téléchargement de Bowtie 0.12.7
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip 

# Décompression de l'archive
RUN unzip /bowtie-0.12.7-linux-x86_64.zip

# Ajout du chemin vers Bowtie à la variable PATH
ENV PATH="$PATH:/bowtie-0.12.7"

# Suppression des fichiers temporaires et inutiles
RUN rm -f /tmp/bowtie-0.12.7-linux-x86_64.zip

