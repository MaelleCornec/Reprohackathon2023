# Spécification de l'image de base Ubuntu
FROM ubuntu:16.04

# Installation des dépendances nécessaires pour DESeq2
RUN apt-get update && apt-get install -y \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Installation de BiocManager et spécification
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/'); \
          BiocManager::install(version='3.2')"

# Installation de DESeq2 
RUN R -e "BiocManager::install('DESeq2')"

# Démarrez R à l'intérieur du conteneur
CMD ["R"]
