# Spécification de l'image Ubuntu comme base
FROM ubuntu:22.04 

# Mise à jour du système et installation des dépendances nécessaires
RUN apt-get update -y && \
    apt-get install -y wget python3 python3-pip

# Téléchargement de Trim Galore version 0.6.4 depuis le site officiel
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz

# Décompression de l'archive
RUN tar xzf 0.6.4.tar.gz

# Ajout du répertoire contenant Trim Galore à la variable PATH
ENV PATH="/TrimGalore-0.6.4:${PATH}"

# Suppresion des fichiers temporaires et inutiles
RUN rm -rf 0.6.4.tar.gz

# Installation de Cutadapt version 1.11
RUN pip3 install cutadapt==1.11

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["trim_galore"]