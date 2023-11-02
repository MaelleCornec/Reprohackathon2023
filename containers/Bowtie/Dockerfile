#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update -y && apt-get install -y \
    build-essential \
    wget

# Téléchargement de Bowtie Version 0.12.7 depuis la source officielle pour l'alignement des séquences
#WORKDIR Reprohackathon2023/containers/Bowtie
#RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-src.zip
#RUN unzip bowtie-0.12.7-src.zip
#RUN rm bowtie-0.12.7-src.zip
#WORKDIR /Reprohackathon2023/containers/Bowtie

# Compilation de Bowtie
#RUN make

# Ajout de Bowtie au PATH
#ENV PATH="/opt/bowtie-0.12.7:${PATH}"

# Définition du répertoire de travail par défaut
#WORKDIR /Reprohackathon2023/containers/Bowtie

# Pour vérifier que Bowtie est installé correctement
#CMD ["bowtie", "--version"]
#CMD ["/bin/bash"]

# Utilisez une image de base Ubuntu comme point de départ
#FROM ubuntu:latest

# Mettez à jour le système et installez les dépendances nécessaires
RUN apt-get update && apt-get install -y gcc make wget build-essential unzip

# Créez un répertoire de travail /app dans l'image
WORKDIR /home/ubuntu/Reprohackathon2023/containers/Bowtie

# Téléchargez Bowtie 0.12.7 depuis le site officiel
RUN wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.7/bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip

# Compilez Bowtie sans changer de répertoire
RUN ls
RUN cd bowtie-0.12.7 && make


# Supprimez les fichiers temporaires et inutiles
#RUN rm -rf bowtie-0.12.7 bowtie-0.12.7-src.zip

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["/home/ubuntu/Reprohackathon2023/containers/Bowtie/bowtie-0.12.7"]


