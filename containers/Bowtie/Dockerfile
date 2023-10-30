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

####### Téléchargement de Bowtie Version 0.12.7 depuis la source officielle pour l'alignement des séquences
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
CMD ["bowtie", "--version"]