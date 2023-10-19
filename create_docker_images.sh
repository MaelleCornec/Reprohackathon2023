#!/bin/bash

# Créez l'image Docker à partir du Dockerfile
docker build -t bowtie:0.12.7 -f Dockerfile .

# Exemple : lancez un conteneur interactif basé sur l'image nouvellement créée
docker run -it bowtie:0.12.7

# Vous pouvez personnaliser davantage ce script en fonction de vos besoins.
