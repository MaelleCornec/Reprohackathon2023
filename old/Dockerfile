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

#############################################################################################
####### SRA Toolkit pour le téléchargement des fichiers fastq #########
# Téléchargement et installation de SRA Toolkit
RUN mkdir docker_SRAToolkit
WORKDIR /docker_SRAToolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    cd sratoolkit.current-ubuntu64 && \
    ./install
# Nettoyage des fichiers temporaires
RUN rm -rf sratoolkit.current-ubuntu64.tar.gz sratoolkit.current-ubuntu64
# Commande par défaut à exécuter lorsque le conteneur démarre
CMD [ "fastq-dump", "--version" ]


#############################################################################################
####### Bowtie Version 0.12.7 pour l'alignement des séquences ######
# Téléchargement de Bowtie 0.12.7 depuis la source officielle
RUN mkdir docker_bowtie
WORKDIR /docker_bowtie
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
RUN unzip download
RUN touch bowtie-build.sh
RUN echo bowtie-build 
#WORKDIR /docker_bowtie/bowtie-0.12.7
#RUN chmod +rwx bowtie-build
#RUN ./bowtie-0.12.7/bowtie-build ../reference.fasta ../reference.gff
# Définissez un répertoire de travail
#WORKDIR ~/Reprohackathon2023/docker_Bowtie_0.12.7
# Par exemple, pour vérifier que Bowtie est installé correctement
#CMD ["bowtie", "--version"]


#############################################################################################
##### CUTADAPT version 1.11 pour le nettoyage des séquences #####
RUN mkdir doker_cutadapt
WORKDIR /docker_cutadapt
RUN pip3 install cutadapt==1.11


#############################################################################################
######### featureCounts version 1.4.6-p3 du package Subreads (paramètres : -t gene -g ID -s 1) 
#####pour le comptage des gènes ######
# Spécification de la version de FeatureCounts désirée (1.4.6-p3)
RUN mkdir docker_featureCounts
WORKDIR /docker_featureCounts
RUN wget http://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    tar -zxvf subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    mv subread-1.4.6-p3-Linux-x86_64/bin/* /usr/local/bin/
# Nettoyage des fichiers temporaires
RUN rm -rf subread-1.4.6-p3-Linux-x86_64.tar.gz subread-1.4.6-p3-Linux-x86_64
# Commande par défaut exécutée lorsque le conteneur démarre
CMD [ "featureCounts", "--version" ]


#############################################################################################
########### R version version 3.4.1 pour l'analyse des données de comptage  #################
# Ajout du référentiel CRAN pour R
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository "deb http://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)/"
# Installation de R version 3.4.1
RUN apt-get update && apt-get install -y r-base=3.4.1-2trusty
# Installation de Bioconductor et le package DESeq2 version 1.16
# pour la normalisation et l'estimation de la dispersion en utilisant les paramètres par défaut
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.4', ask = FALSE)"
RUN R -e "BiocManager::install('DESeq2', version = '1.16')"
# Nettoyage du cache de R
RUN R -e "q()"