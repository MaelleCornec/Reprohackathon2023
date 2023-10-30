#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update && apt-get install -y python3 python3-pip
RUN apt-get update && apt-get install -y build-essential zlib1g-dev curl unzip
RUN apt install wget
RUN apt-get update && apt-get install -y subread
RUN apt-get update && apt-get install -y software-properties-common
RUN apt-get update && apt-get install -y libxml2 libssl-dev libcurl4-openssl-dev
RUN apt-get install -y wget

# Téléchargement et installation de SRA Toolkit pour le téléchargement des fichiers fastq 

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    cd sratoolkit.current-ubuntu64 && \
    ./install
# Nettoyage des fichiers temporaires
RUN rm -rf sratoolkit.current-ubuntu64.tar.gz sratoolkit.current-ubuntu64
# Commande par défaut à exécuter lorsque le conteneur démarre
CMD [ "fastq-dump", "--version" ]