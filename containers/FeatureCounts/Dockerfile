#Spécification de l'image de base qu'on souhaite utiliser 
FROM ubuntu:22.04 

#Installation des dépendances et des logiciels
RUN apt-get update && apt-get install -y build-essential

# featureCounts version 1.4.6-p3 du package Subreads (paramètres : -t gene -g ID -s 1) pour le comptage des gènes ######
RUN wget http://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    tar -zxvf subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    mv subread-1.4.6-p3-Linux-x86_64/bin/* /usr/local/bin/
# Nettoyage des fichiers temporaires
RUN rm -rf subread-1.4.6-p3-Linux-x86_64.tar.gz subread-1.4.6-p3-Linux-x86_64

# Commande par défaut exécutée lorsque le conteneur démarre
CMD [ "featureCounts", "--version" ]

