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

# R version version 3.4.1 pour l'analyse des données de comptage 
# Ajout du référentiel CRAN pour R
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository "deb http://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)/"
# Installation de R version 3.4.1
RUN apt-get update && apt-get install -y r-base=3.4.1-2trusty
# Installation de Bioconductor et le package DESeq2 version 1.16 pour la normalisation et l'estimation de la dispersion en utilisant les paramètres par défaut
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.4', ask = FALSE)"
RUN R -e "BiocManager::install('DESeq2', version = '1.16')"
# Nettoyage du cache de R
RUN R -e "q()"