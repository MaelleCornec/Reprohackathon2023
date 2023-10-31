#!/bin/bash

# On donne tous les droits aux utilisateurs 
chmod -R 777 /home/ubuntu/Reprohackathon2023/containers

# Création de l'image sratoolkit à partir du Dockerfile dans le répertoire containers/SRA_Toolkit
/home/ubuntu/Reprohackathon2023/containers/SRA_Toolkit/Dockerfile
docker build -t sratoolkit -f Dockerfile .
# Exécution du conteneur sratoolkit interactif 
docker run -it sratoolkit

# Création de l'image trimgalore à partir du Dockerfile dans le répertoire containers/TrimGalore
/home/ubuntu/Reprohackathon2023/containers/TrimGalore/Dockerfile
docker build -t trimgalore -f Dockerfile .
# Exécution du conteneur trimgalore interactif 
docker run -it trimgalore

# Création de l'image bowtie à partir du Dockerfile dans le répertoire containers/Bowtie
/home/ubuntu/Reprohackathon2023/containers/Bowtie/Dockerfile
docker build -t bowtie:0.12.7 -f Dockerfile .
# Exécution du conteneur bowtie interactif 
docker run -it bowtie:0.12.7

# Création de l'image featurecounts à partir du Dockerfile dans le répertoire containers/FeatureCounts
/home/ubuntu/Reprohackathon2023/containers/FeatureCounts/Dockerfile
docker build -t featurecounts -f Dockerfile .
# Exécution du conteneur featurecounts interactif 
docker run -it featurecounts

# Création de l'image r:3.4.1 à partir du Dockerfile dans le répertoire containers/R
/home/ubuntu/Reprohackathon2023/containers/R/Dockerfile
docker build -t r:3.4.1 -f Dockerfile .
# Exécution du conteneur r interactif 
docker run -it r:3.4.1
