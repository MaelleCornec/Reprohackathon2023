# Utilisez une image de base qui prend en charge Bowtie
FROM ubuntu:latest

# Mettez à jour le système de base et installez les dépendances nécessaires, y compris curl
RUN apt-get update && apt-get install -y build-essential zlib1g-dev curl unzip
RUN apt install wget

# Téléchargez Bowtie 0.12.7 depuis la source officielle et compilez-le
RUN mkdir docker_bowtie
WORKDIR /docker_bowtie

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download
RUN unzip download

RUN touch bowtie-build.sh
RUN echo bowtie-build 

#WORKDIR /docker_bowtie/bowtie-0.12.7
#RUN chmod ugo+rwx bowtie-build

#RUN ./bowtie-0.12.7/bowtie-build ../reference.fasta ../reference.gff

# Définissez un répertoire de travail
#WORKDIR ~/Reprohackathon2023/docker_Bowtie_0.12.7

# Par exemple, pour vérifier que Bowtie est installé correctement
#CMD ["bowtie", "--version"]
