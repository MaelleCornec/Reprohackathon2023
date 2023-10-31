# Spécification de l'image Ubuntu 22.04 comme base
FROM ubuntu:22.04

# Installation des dépendances nécessaires
RUN apt-get update && apt-get install -y wget

# Téléchargement de SRA Toolkit depuis NCBI
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.2.10.8-ubuntu64.tar.gz && \
    mv sratoolkit.2.10.8-ubuntu64 /opt/sratoolkit && \
    rm sratoolkit.2.10.8-ubuntu64.tar.gz

# Ajout du répertoire de SRA Toolkit à la variable PATH
ENV PATH="/opt/sratoolkit/bin:${PATH}"

# Définition du répertoire de travail
WORKDIR /home/ubuntu/Reprohackathon2023/data

# Commande par défaut (exécutez `vdb-config -version` pour vérifier l'installation)
CMD ["vdb-config", "-version"]