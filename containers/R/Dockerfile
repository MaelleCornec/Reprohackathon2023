#Spécification de l'image de base qu'on souhaite utiliser 
#FROM ubuntu:22.04 

# R version version 3.4.1 pour l'analyse des données de comptage 
# Ajout du référentiel CRAN pour R
#RUN apt-get update && apt-get install -y software-properties-common
#RUN add-apt-repository "deb http://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)/"
# Installation de R version 3.4.1
#RUN apt-get update && apt-get install -y r-base=3.4.1-2trusty
# Installation de Bioconductor et le package DESeq2 version 1.16 pour la normalisation et l'estimation de la dispersion en utilisant les paramètres par défaut
#RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
#RUN R -e "BiocManager::install(version = '3.4', ask = FALSE)"
#RUN R -e "BiocManager::install('DESeq2', version = '1.16')"
# Nettoyage du cache de R
#RUN R -e "q()"



# Mettez à jour le système et installez les dépendances nécessaires
#RUN apt-get update && apt-get install -y wget

# Téléchargez R version 3.4.1 depuis un miroir CRAN
#RUN wget https://cran.r-project.org/src/base/R-3/R-3.4.1.tar.gz

# Décompressez l'archive
#RUN tar xzf R-3.4.1.tar.gz

# Compilez R depuis les sources
#RUN cd R-3.4.1 && ./configure && make && make install

# Supprimez les fichiers temporaires et inutiles
#RUN rm -rf R-3.4.1 R-3.4.1.tar.gz

# Installez des dépendances pour Bioconductor et DESeq2
#RUN apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev libcairo2-dev

# Installez Bioconductor
#RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org/'); BiocManager::install()"

# Installez le package DESeq2 version 1.16 depuis Bioconductor
#RUN R -e "BiocManager::install('DESeq2', version = '1.16')"

# Commande par défaut à exécuter lorsque le conteneur démarre
#CMD ["R"]

# Utilisez une image de base basée sur R
FROM rocker/r-ver:3.4.1

# Installez les dépendances nécessaires pour Bioconductor et DESeq2
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libcairo2-dev

# Installez Bioconductor
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org/'); BiocManager::install()"

# Installez le package DESeq2 version 1.16 depuis Bioconductor
RUN R -e "BiocManager::install('DESeq2', version = '1.16')"

# Commande par défaut à exécuter lorsque le conteneur démarre
CMD ["R"]
