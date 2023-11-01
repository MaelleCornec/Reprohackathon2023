#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update -y && apt-get install -y \
    build-essential \
    wget

# Téléchargement de Bowtie Version 0.12.7 depuis la source officielle pour l'alignement des séquences
WORKDIR /opt
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip
RUN rm bowtie-0.12.7-src.zip
WORKDIR /opt/bowtie-0.12.7

# Compilation de Bowtie
RUN make

# Ajout de Bowtie au PATH
ENV PATH="/opt/bowtie-0.12.7:${PATH}"

# Définition du répertoire de travail par défaut
WORKDIR /opt/bowtie-0.12.7

# Pour vérifier que Bowtie est installé correctement
CMD ["bowtie", "--version"]
#CMD ["/bin/bash"]