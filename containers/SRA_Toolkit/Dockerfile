#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update && apt-get install -y wget

# Téléchargement et installation de SRA Toolkit pour le téléchargement des fichiers fastq 
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.2.10.8-ubuntu64.tar.gz && \
    mv sratoolkit.2.10.8-ubuntu64 /opt/sratoolkit && \
    rm sratoolkit.2.10.8-ubuntu64.tar.gz

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD [ "fastq-dump", "--version" ]