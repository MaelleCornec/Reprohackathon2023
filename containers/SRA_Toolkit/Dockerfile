#Spécification de l'image de base 
FROM ubuntu:18.04

#Téléchargament des dépendances nécessaires 
RUN apt-get update \
	&& apt-get install -y wget libxml-libxml-perl

#Téléchargement de SRA Toolkit version 2.10.0 depuis la source officielle 
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz -O /tmp/sratoolkit.tar.gz \
	&& tar zxvf /tmp/sratoolkit.tar.gz -C /opt/ && rm /tmp/sratoolkit.tar.gz

ENV PATH="/opt/sratoolkit.2.10.0-ubuntu64/bin/:${PATH}"